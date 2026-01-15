import re
import os
from collections import Counter
from pathlib import Path

# ============================ VULCAN NETWORK ==============================================

def vulcan(path, allowed_elements, metadata="N"):
    print("FILTERING VULCAN CHEMICAL NETWORK....")

    def parse_formula(formula):
        tokens = re.findall(r'([A-Z][a-z]?)(\d*)', formula)
        comp = Counter()
        for el, n in tokens:
            comp[el.upper()] += int(n) if n else 1
        return comp

    def species_allowed(species):
        if species.strip().upper() == "M":
            return True
        species_clean = re.sub(r'_\d+', '', species.strip())
        species_atoms = parse_formula(species_clean)
        allowed_upper = {el.upper() for el in allowed_elements}
        return set(species_atoms.keys()).issubset(allowed_upper)

    def extract_species(reaction_str):
        try:
            lhs, rhs = reaction_str.split("->")
        except ValueError:
            return []
        parts = lhs.split("+") + rhs.split("+")
        return [re.sub(r'_\d+', '', p.strip()) for p in parts]

    filename = "VULCAN_output_" + "_".join(sorted(el for el in allowed_elements)) + ".dat"
    kept = []

    with open(path, "r") as f:
        for line in f:
            raw = line.rstrip()
            if not raw:
                continue

            if raw.startswith("#"):
                if metadata.upper() == "Y":
                    kept.append(raw)
                continue

            m = re.search(r'\[(.*?)\]', raw)
            if not m:
                continue

            reaction_str = m.group(1)
            species = extract_species(reaction_str)

            if species and all(species_allowed(s) for s in species):
                kept.append(raw)

    if not kept:
        print("No reactions matched the allowed elements.")
        return

    with open(filename, "w") as out:
        for line in kept:
            out.write(line + "\n")

    print(f"Filtered reactions written to {filename}\n")




# ============================ HU.  NETWORK ================================================

def hu(path, allowed_elements, metadata="N"):
    print("FILTERING HU CHEMICAL NETWORK....")

    # ---------- helpers ----------

    def normalize_species(species):
        """
        Remove excited-state notation and labels.
        Examples:
        O(^1D), O^1D, O(1D) -> O
        O_2 -> O2
        """
        s = species.strip()

        # remove everything from first ^ onward
        s = re.split(r'\^', s)[0]

        # remove parenthetical states
        s = re.sub(r'\(.*?\)', '', s)

        # remove underscore labels like O_2
        s = re.sub(r'_\d+', '', s)

        return s.strip()

    def parse_formula(formula):
        """Return atomic composition as uppercase keys"""
        tokens = re.findall(r'([A-Z][a-z]?)(\d*)', formula)
        comp = Counter()
        for el, n in tokens:
            comp[el.upper()] += int(n) if n else 1
        return comp

    def species_allowed(species, allowed_elements):
        """
        Check if species is made only of allowed elements.
        Excited states already collapsed.
        """
        if species.strip().upper() == "M":
            return True

        clean = normalize_species(species)
        if not clean:
            return False

        atoms = parse_formula(clean)
        allowed_upper = {el.upper() for el in allowed_elements}

        return set(atoms.keys()).issubset(allowed_upper)

    def extract_hu_species(line):
        """
        Extract species from Hu-format reaction line.
        Example:
        R3   C + H2         CH + H
        """
        text = line[4:36]
        if not text.strip:
            return []

        text = line[4:36]
        if not text.strip():
            print("no text here, empty line")
            return []

        parts = re.split(r'\s{2,}', text, maxsplit=1)
        if len(parts) != 2:
            print("Less than two reactants and products here")
            return[]

        reactants, products = parts

        species = []

        for side in (reactants, products):
            species.extend(s.strip() for s in side.split("+"))
        
        return [s for s in species if s]

    def extract_photo_species(line):
        """
        Extract species from photolysis reaction first line:
        1  O_2  O + O   ...
        """
        parts = line.split()
        if len(parts) < 4:
            return []

        reactant = parts[1]
        products = parts[2].split("+")

        return [reactant] + [p.strip() for p in products]

    # ---------- main logic ----------

    tag = "_".join(sorted(el for el in allowed_elements))
    outfile = f"Hu_output_{tag}.dat"

    kept = []
    mode = "hu"
    current_row = None
    photo_keep_current = False

    with open(path, "r") as f:
        for line in f:
            raw = line.rstrip("\n")

            # metadata
            if raw.startswith("#"):
                if metadata.upper() == "Y":
                    kept.append(raw)
                if raw.strip().startswith("#Table 2"):
                    mode = "photo"
                continue

            if not raw.strip():
                continue

            # HU reactions
            if mode == "hu":
                if not raw.lstrip()[0] in ("R", "M", "T"):
                    continue

                species = extract_hu_species(raw)
                if species and all(species_allowed(s, allowed_elements) for s in species):
                    kept.append(raw)

            # Photolysis reactions
            elif mode == "photo":
                parts = raw.split()

                if parts and parts[0].isdigit():
                    if current_row and photo_keep_current:
                        kept.append(current_row)

                    species = extract_photo_species(raw)
                    photo_keep_current = (
                        species and
                        all(species_allowed(s, allowed_elements) for s in species)
                    )
                    current_row = raw
                else:
                    if current_row and photo_keep_current:
                        current_row += "\n" + raw

        if current_row and photo_keep_current:
            kept.append(current_row)

    if not kept:
        print("No reactions matched the allowed elements.")
        return

    with open(outfile, "w") as out:
        for block in kept:
            out.write(block + "\n")

    print(f"Filtered reactions written to {outfile}\n")




#=================================== MOSES NETWORK ==========================================

def moses(path, allowed_elements, metadata="N", info = "N"):
    """
    Master function for Moses chemical network.
    path: folder containing ReadMe and TableS2 .txt files
    allowed_elements: list of allowed elements, e.g., ['O', 'H', 'C']
    metadata: 'Y' to include ReadMe, else 'N'
    """

    print("FILTERING MOSES CHEMICAL NETWORK....")

    # ---------------- Helper Functions ----------------
    def parse_formula(formula):
        """Return atomic composition as uppercase keys"""
        tokens = re.findall(r'([A-Z][a-z]?)(\d*)', formula)
        comp = Counter()
        for el, n in tokens:
            comp[el.upper()] += int(n) if n else 1
        return comp

    def species_allowed(species, allowed_elements):
        """
        Check species composition against allowed_elements
        Excited states stripped
        'M' always allowed
        """
        species = strip_excited_state(species)

        if species.upper() == "M":
            return True

        species_atoms = parse_formula(species)
        allowed_upper = {el.upper() for el in allowed_elements}
        return set(species_atoms.keys()).issubset(allowed_upper)

    def strip_excited_state(species):
        """
        Collapse excited states:
        O(1D) -> O
        (1)CH2, (3)CH2 -> CH2
        """
        # Remove parenthesis with digits for CH2
        species = re.sub(r'\(\d+\)', '', species)
        # Remove O(1D) style
        species = re.sub(r'O\(\d+D\)', '', species)
        
        return species.strip()

    def extract_species(line):
        if '=' not in line:
            return []

        # remove reaction number
        line = re.sub(r'^\s*\d+\)\s*', '', line)

        # cut off kinetics (everything starting with k, k0, koo)
        line = re.split(r'\bk0\b|\bkoo\b|\bk\b', line)[0]

        if '=' not in line:
            return[]

        lhs, rhs = line.split('=', 1)

        parts = lhs.split('+') + rhs.split('+')
        return [p.strip() for p in parts if p.strip()]


    # ---------------- File Setup ----------------
    tag = "_".join(sorted(el for el in allowed_elements))
    parent_folder = os.path.dirname(path)
    output_file = os.path.join(parent_folder, 
                               f"MOSES_output_{tag}.dat")
    kept = []

    # ---------------- Metadata ----------------
    readme_file = os.path.join(path, "ReadMe")
    if metadata.upper() == "Y" and os.path.exists(readme_file):
        with open(readme_file, "r") as f:
            for line in f:
                kept.append(f"# {line.rstrip()}")
        kept.append("\n\n\n\n")

    # Add fixed comments
    if info == 'Y':
        kept.append("## IMPORTANT INFORMATION (added by Junaid)")
        kept.append("# For rate coefficients please check TableS1.pdf")
        kept.append("# R1-137 are photolysis reactions")

    # ---------------- Read TableS2 Text File ----------------
    txt_files = [f for f in os.listdir(path) if f.endswith(".txt")]
    if not txt_files:
        print("No .txt files found in folder.")
        return

    # Use the specified file
    txt_file = os.path.join(path, "TableS2_hd189_1d2pi.txt")
    if not os.path.exists(txt_file):
        txt_file = os.path.join(path, txt_files[0])  # fallback
        print(f"Specified file not found, using {txt_file}")

    # ---------------- Process File ----------------
    stop_line_pattern = "CONCENTRATIONS OF"
    reactions_started = False

    with open(txt_file, "r") as f:
        for lineno, line in enumerate(f, 1):
            raw = line.rstrip("\n")
            
            # First 29 characters as header
            if info == 'Y':
                if lineno == 1:
                    kept.append("\n\n\n")
                    kept.append(raw[:29])
                    continue


            # Stop condition
            if stop_line_pattern in raw or lineno > 2115:
                break

            # Skip until "REACTIONS:" is found
            if not reactions_started:
                if info == 'Y':
                    kept.append(raw)
                if "REACTIONS:" in raw:
                    reactions_started = True
                continue

            # Extract reaction number
            if re.match(r'^\s*\d+\)', raw):
                species = extract_species(raw)
                if species and all(species_allowed(s, allowed_elements) for s in species):
                    kept.append(raw)
            else:
                species = extract_species(raw)
                if species and all(species_allowed(s, allowed_elements) for s in species):
                    kept[-1] += "\n" + raw

    # ---------------- Write Output ----------------
    if not kept:
        print("No reactions matched the allowed elements.")
        return

    with open(output_file, "w") as out:
        for block in kept:
            out.write(block + "\n")

    print(f"Filtered reactions written to {output_file}\n")




#=================================== STAND NETWORK ==========================================

def stand(path, allowed_elements, metadata="N"):
    """
    Master function for STAND chemical network.

    Parameters
    ----------
    path : str
        Path to the STAND input file
    allowed_elements : list
        Elements to keep, e.g. ['H', 'C', 'O']
    metadata : str
        'Y' to include metadata, 'N' otherwise
    """

    print("FILTERING STAND2020 CHEMICAL NETWORK....")

    allowed_upper = {el.upper() for el in allowed_elements}

    # -------------------------------------------------
    # Helpers
    # -------------------------------------------------
    SPECIES_ALIASES = {
                        "Oxyrane": "C2H4O"
                    }

    def parse_formula(formula):
        """
        Parse a plain chemical formula into elements.
        Example: CH3O2 -> {'C':1,'H':3,'O':2}
        """
        tokens = re.findall(r'([A-Z][a-z]?)(\d*)', formula)
        comp = Counter()
        for el, n in tokens:
            comp[el.upper()] += int(n) if n else 1
        return comp

    def reduced_formula(species):
        """
        Reduce a species ONLY for element checking.
        Output is never written to file.
        """

        s = species.strip()

        # Always-allowed helpers
        if re.fullmatch(r'M', s, re.IGNORECASE):
            return None
        if re.search(r'gamma', s, re.IGNORECASE):
            return None
        if re.search(r'e', s) and '-' in s:
            return None
         
        if s in SPECIES_ALIASES:
            return SPECIES_ALIASES[s]

        # Remove molecular state: O(^1D), O_2 (^1Δ), etc
        s = re.sub(r'\([^)]*\)', '', s)

        # Remove charge / excitation notation
        s = re.sub(r'[\^\+\-\*]', '', s)

        # Remove underscores only for parsing
        s = s.replace('_', '')

        return s.strip()

    def species_allowed(species):
        reduced = reduced_formula(species)

        # M, gamma, electrons
        if reduced is None:
            return True

        atoms = parse_formula(reduced)
        return set(atoms.keys()).issubset(allowed_upper)

    def extract_species(line):
        """
        Extract species from a STAND reaction line.
        Reaction separator is => .
        Kinetics are ignored automatically.
        """

        if "=>" not in line:
            return []

        # Remove anything after kinetics start (units, numbers, refs)
        reaction_part = re.split(r'\s{2,}', line)[0]

        try:
            lhs, rhs = reaction_part.split("=>", 1)
        except ValueError:
            return []

        parts = lhs.split("+") + rhs.split("+")
        return [p.strip() for p in parts if p.strip()]

    
    # -------------------------------------------------
    # Output setup
    # -------------------------------------------------

    tag = "_".join(sorted(allowed_upper))
    out_file = os.path.join(
        os.path.dirname(path),
        f"stand_output_{tag}.dat"
    )

    kept = []

    # -------------------------------------------------
    # Read file
    # -------------------------------------------------

    metadata_done = False
    metadata_end_pattern = r'Krasnopolsky \(2007Icar..191...25K\); Zhang et al. \(2012Icar..217..714Z\).'
    

    with open(path, "r") as f:
        kept.append(("#      REACTION                                      ALPHA" +        
          "         BETA    GAMMA  REF     (added by JRB)"))
        for lineno, line in enumerate(f, 1):
            raw = line.rstrip("\n")

            # Metadata block
            if not metadata_done:
                if re.search(metadata_end_pattern, raw): #or "=>" in raw:
                    metadata_done = True
                    if metadata.upper() == "Y":
                        kept.append(raw)
                    continue
                

                if metadata.upper() == "Y":
                    kept.append(raw)
                continue
                

            # Empty lines
            if not raw.strip():
                continue

            # Extract and filter reactions
            species = extract_species(raw)

            if species and all(species_allowed(s) for s in species):
                kept.append(raw)
        
        
    # -------------------------------------------------
    # Write output
    # -------------------------------------------------

    if not kept:
        print("No reactions matched the allowed elements.")
        return

    with open(out_file, "w") as out:
        for line in kept:
            out.write(line + "\n")

    print(f"Filtered reactions written to {out_file}\n")




#=================================== VELLIET_VENOT NETWORK ==========================================

def velliet_venot(path, allowed_elements, metadata="N"):
    """
    Master function for the velliet_venot chemical network.
    """

    print("FILTERING VELLIET_VENOT CHEMICAL NETWORK...")

    allowed_upper = {el.upper() for el in allowed_elements}

    # -------------------------------------------------
    # Helpers
    # -------------------------------------------------

    def parse_formula(formula):
        """
        Parse chemical formula into elements.
        Example: C2H4O -> {'C':2,'H':4,'O':1}
        """
        tokens = re.findall(r'([A-Z][a-z]?)(\d*)', formula)
        comp = Counter()
        for el, n in tokens:
            comp[el.upper()] += int(n) if n else 1
        return comp

    def reduce_species(species):
        """
        Reduce species name ONLY for element checking.
        Never written to output.
        """

        s = species.strip()

        # Always allowed
        if s == "HV":
            return None

        # Remove brackets: CH2CHN(S)
        s = re.sub(r'\([^)]*\)', '', s)

        # Remove leading numeric state: 1CH2, 3CH2
        s = re.sub(r'^\d+', '', s)

        # Remove vibrational / electronic state letters
        # O1D, O3P, OHV, N2D, S1D
        s = re.sub(r'(D|P|V)$', '', s)

        # Remove hyphen states: C3H4-A, C2O5H5-1
        s = re.sub(r'-.*$', '', s)

        # Remove cyclic prefixes: c-C2H3N
        s = re.sub(r'^c-', '', s)

        return s.strip()

    def species_allowed(species):
        reduced = reduce_species(species)

        # HV or empty
        if reduced is None or not reduced:
            return True

        atoms = parse_formula(reduced)
        return set(atoms.keys()).issubset(allowed_upper)

    def extract_species_from_line(line):
        """
        velliet_venot format:
        reactants are in fixed-width columns before products.
        We split on large whitespace blocks.
        """
        parts = re.split(r'\s{2,}', line.strip())

        # Heuristic: first half are reactants, next are products
        # This works consistently for VV files
        species = []
        for p in parts:
            if re.search(r'[A-Z]', p):
                species.append(p)
        return species

    # -------------------------------------------------
    # Output setup
    # -------------------------------------------------

    tag = "_".join(sorted(allowed_upper))
    out_file = os.path.join(os.path.dirname(path), f"velliet_venot_output_{tag}.dat")

    kept = []

    # -------------------------------------------------
    # Read README
    # -------------------------------------------------

    readme = os.path.join(path, "readme.txt")
    if metadata.upper() == "Y" and os.path.exists(readme):
        kept.append("readme.txt\n")
        with open(readme, "r") as f:
            for line in f:
                kept.append(line.rstrip("\n"))
        kept.append("")

    # -------------------------------------------------
    # File order
    # -------------------------------------------------

    files = ["photodissociations.dat"] + [
        f"reactions_{i}.dat" for i in range(1, 14)
    ]

    # -------------------------------------------------
    # Process files
    # -------------------------------------------------

    for fname in files:
        fpath = os.path.join(path, fname)
        if not os.path.exists(fpath):
            continue

        file_hits = []

        with open(fpath, "r") as f:
            for line in f:
                raw = line.rstrip("\n")

                if not raw.strip():
                    continue

                species = extract_species_from_line(raw)

                if species and all(species_allowed(s) for s in species):
                    file_hits.append(raw)

        if file_hits:
            kept.append(fname)
            kept.append("-" * len(fname))
            kept.extend(file_hits)
            kept.append("")

    # -------------------------------------------------
    # Write output
    # -------------------------------------------------

    if not kept:
        print("No reactions matched the allowed elements.")
        return

    with open(out_file, "w") as out:
        for line in kept:
            out.write(line + "\n")

    print(f"Filtered reactions written to {out_file}\n")



#=================================== AGUNDEZ NETWORK ==========================================


def agundez(path, allowed_elements, metadata="N"):
    print("FILTERING AGUNDEZ CHEMICAL NETWORK....")

    allowed_upper = {el.upper() for el in allowed_elements}

    # -------------------------------------------------
    # Helpers
    # -------------------------------------------------
    def parse_formula(formula):
        tokens = re.findall(r'([A-Z][a-z]?)(\d*)', formula)
        comp = Counter()
        for el, n in tokens:
            comp[el.upper()] += int(n) if n else 1
        return comp

    def reduced_formula(species):
        s = species.strip()
        if s.upper() == "M":
            return None
        s = s.split("_")[0]  # excited states
        s = re.sub(r'[^A-Za-z0-9]', '', s)
        return s.strip()

    def species_allowed(species):
        reduced = reduced_formula(species)
        if reduced is None or not reduced:
            return True
        atoms = parse_formula(reduced)
        return set(atoms.keys()).issubset(allowed_upper)

    def extract_species(line):
        if "=" not in line:
            return []
        reaction_part = line.split(":", 1)[0]
        try:
            lhs, rhs = reaction_part.split("=", 1)
        except ValueError:
            return []
        parts = lhs.split("+") + rhs.split("+")
        return [p.strip() for p in parts if p.strip()]

    # -------------------------------------------------
    # Output setup
    # -------------------------------------------------
    tag = "_".join(sorted(allowed_upper))
    out_file = os.path.join(
        os.path.dirname(path),
        f"Agundez_output_{tag}.dat"
    )

    kept = []
    discard_kept = []

    # -------------------------------------------------
    # State flags
    # -------------------------------------------------
    in_metadata = True
    in_bimolecular_block = False
    discard_block = False

    # -------------------------------------------------
    # Main loop
    # -------------------------------------------------
    with open(path, "r") as f:
        for line in f:
            raw = line.rstrip("\n")

            # ----- Block markers -----
            if raw.startswith("! Bimolecular reactions"):
                in_metadata = False
                in_bimolecular_block = True
                if metadata.upper() == "Y":
                    kept.append(raw)
                continue

            if raw.startswith("! Reactions discarded because they involve species not included"):
                in_bimolecular_block = False
                discard_block = True
                discard_kept = []
                if metadata.upper() == "Y":
                    kept.append(raw)
                continue

            if raw.startswith("! # References"):
                in_metadata = True
                discard_block = False
                if metadata.upper() == "Y":
                    kept.append(raw)
                continue

            # ----- Inside BIMOL block -----
            if in_bimolecular_block and raw.startswith("!"):
                # "!<space>" → metadata
                if re.match(r'!\s', raw):
                    if metadata.upper() == "Y":
                        kept.append(raw)
                    continue

                # "!<no space>" → reaction
                if re.match(r'!\S', raw):
                    reaction_line = raw[1:].lstrip()
                    species = extract_species(reaction_line)
                    if species and all(species_allowed(s) for s in species):
                        kept.append(raw)
                    continue

            # ----- Discarded reactions block -----
            if discard_block and "=" in raw:
                species = extract_species(raw)
                if species and all(species_allowed(s) for s in species):
                    discard_kept.append(raw)
                continue

            # ----- Normal metadata -----
            if raw.startswith("!") and in_metadata:
                if metadata.upper() == "Y":
                    kept.append(raw)
                continue

            # ----- Normal reactions -----
            if "=" in raw:
                species = extract_species(raw)
                if species and all(species_allowed(s) for s in species):
                    kept.append(raw)

    if discard_kept:
        kept.append("\n!Reactions discarded because they involve species not included\n")
        kept.extend(discard_kept)
    
    if not kept:
        print("No reactions matched the allowed elements.")
        return

    with open(out_file, "w") as out:
        for line in kept:
            out.write(line + "\n")

    print(f"Filtered reactions written to {out_file}")


#######################################################################################################
#######################################################################################################
#######################################################################################################

#================================= FILE REPROCESSING ==================================================


"""
This code is used to identify the unique (many reactions are repeated in different networks) 
reactions in the output files generated using above functions. 
USED IN: chem_sorting.py
"""



def reprocess_velliet_file(path):
    path = Path(path)
    out_path = path.with_name("reprocess_" + path.name)

    processed = []

    with open(path) as f:
        for line in f:
            if not line.strip():
                continue

            if line[0].islower() or line.startswith("-"):
                continue

            line = line[:110].rstrip()

            # split into text and space blocks
            parts = re.split(r"(\s+)", line)

            new = []
            for p in parts:
                if p.isspace():
                    n = len(p)
                    if n>=3 and n <= 11:
                        new.append(" + ")
                    elif n >= 30:
                        new.append(" = ")
                    else:
                        new.append(" ")
                else:
                    new.append(p.strip())

            reaction = "".join(new)

            # clean up accidental doubles
            reaction = re.sub(r"\s+\+\s+\+", " + ", reaction)
            reaction = re.sub(r"\s+=\s+=", " = ", reaction)
            reaction = re.sub(r"\s+", " ", reaction)

            processed.append(reaction.strip())

    with open(out_path, "w") as f:
        for r in processed:
            f.write(r + "\n")
    print(f"Done reprocessing {path}\n")
    return out_path


def reprocess_agundez_file(input_path):
    """
    Reprocess an Agundez_output_*.dat file:
      - keep only first 50 characters of each line
      - convert tokens like A_BC -> A(BC), where BC goes until next whitespace
      - output is exactly 50 characters wide per line (padded or trimmed)
      - writes to reprocess_<original_name>.dat
    """
    input_path = Path(input_path)
    if input_path.suffix.lower() != ".dat":
        raise ValueError("Expected a .dat input file")

    out_path = input_path.with_name("reprocess_" + input_path.name)

    # Match a "word" chunk containing underscore, stopping at whitespace.
    # Example matches: A_BC, A_3X, O_1D, CH3_foo
    pattern = re.compile(r"([^\s_]+)_([^\s]+)")

    with open(input_path, "r") as fin, open(out_path, "w") as fout:
        for line in fin:
            # Take only the first 50 columns from the original line
            chunk = line[:50].rstrip("\n")

            # Replace underscore form within that 50-col chunk
            def repl(m):
                left = m.group(1)
                right = m.group(2)
                return f"{left}({right})"

            chunk = pattern.sub(repl, chunk)

            # Enforce exact width 50 after replacement
            if len(chunk) < 50:
                chunk = chunk.ljust(50)
            else:
                chunk = chunk[:50]

            fout.write(chunk + "\n")
    print(f"Done reprocessing {input_path}\n")
    return str(out_path)


def reprocess_hu_file(input_path):

    input_path = Path(input_path)
    if input_path.suffix.lower() != ".dat":
        raise ValueError("Expected a .dat input file")

    out_path = input_path.with_name("reprocess_" + input_path.name)

    first_gap_re = re.compile(r"\s{2,}")

    def normalize_caret_notation_in_line(line: str) -> str:
        tokens = line.split()
        out = []

        for tok in tokens:
            if "^" not in tok:
                out.append(tok)
                continue

            # Base is everything before the first caret, but Hu often has '(' there: O(^1^D)
            base = tok.split("^", 1)[0].rstrip("(").rstrip(")")

            # Pull alphanumerics after carets: ^1 ^D -> ["1","D"] -> "1D"
            exc = re.findall(r"\^([A-Za-z0-9]+)", tok)

            if exc:
                out.append(f"{base}({''.join(exc)})")
            else:
                # fallback: remove carets and clean stray parentheses
                out.append(tok.replace("^", "").replace("((", "(").replace("))", ")"))

        return " ".join(out)

    with open(input_path, "r") as fin, open(out_path, "w") as fout:
        for line in fin:
            if not line.startswith(("R", "M", "T")):
                continue

            starts_with_M = line.startswith("M")

            chunk = line[4:37].rstrip("\n")

            m = first_gap_re.search(chunk)
            if not m:
                continue

            lhs = chunk[:m.start()].rstrip()
            rhs = chunk[m.end():].lstrip()

            if starts_with_M:
                lhs = lhs.strip()
                rhs = rhs.strip()
                if lhs:
                    lhs += " + M"
                if rhs:
                    rhs += " + M"

            out_line = f"{lhs} = {rhs}".strip()
            out_line = normalize_caret_notation_in_line(out_line)

            fout.write(out_line + "\n")
    print(f"Done reprocessing {input_path}\n")
    return str(out_path)


def reprocess_vulcan_file(input_path):
    """
    VULCAN reprocess:

    - For each line, take only cols 7–42 (1-based) => Python slice [6:42]
    - Replace '->' with '='
    - Species rules (token-bounded by whitespace):
        * 'O_1' -> 'O(1D)' (VULCAN only)
        * Any token containing underscores: base_suffix1_suffix2... -> base(suffix1suffix2...)
          Examples:
            - H2O_l_s -> H2O(ls)
            - CH2_1   -> CH2(1)
            - A_a_b   -> A(ab)
            - A_B     -> A(B)
    - Output file name: insert 'reprocess_' before the original filename
      e.g. VULCAN_output_H_He_O.dat -> reprocess_VULCAN_output_H_He_O.dat

    Returns output path as a string.
    """
    input_path = Path(input_path)
    if input_path.suffix.lower() != ".dat":
        raise ValueError("Expected a .dat input file")

    out_path = input_path.with_name("reprocess_" + input_path.name)

    # convert token with underscores to base(suffixes_concatenated)
    def underscore_to_parens(tok: str) -> str:
        tok = tok.strip()
        if tok == "O_1":
            return "O(1D)"
        if "_" not in tok:
            return tok

        base, rest = tok.split("_", 1)

        # For any further underscores, concatenate the pieces (a_b -> ab)
        rest_compact = rest.replace("_", "")

        # If rest is empty (rare), return base
        if not rest_compact:
            return base

        return f"{base}({rest_compact})"

    with open(input_path, "r") as fin, open(out_path, "w") as fout:
        for line in fin:
            # Take cols 7–42 (1-based) => indices 6..41
            chunk = line[6:42].rstrip("\n")

            if "->" not in chunk:
                continue

            # Replace arrow with equals for display
            chunk = chunk.replace("->", "=")

            # Tokenize by whitespace, transform species tokens, then re-join with single spaces
            tokens = chunk.split()
            tokens = [underscore_to_parens(t) for t in tokens]

            # Clean spacing around '='
            out_line = " ".join(tokens)
            out_line = re.sub(r"\s*=\s*", " = ", out_line).strip()

            fout.write(out_line + "\n")
    print(f"Done reprocessing {input_path}\n")
    return str(out_path)


def reprocess_moses_file(input_path):
    """
    Moses reprocess:

    - Read columns 7–70 (1-based) => Python slice [6:70]
    - Normalize spacing around '+' and '='
    - Expand stoichiometric prefixes like 2H, 2OH, 2HO2, 3O, etc.
      ONLY when the number is a PREFIX and the species starts with a capital letter.
      Molecules like H2O are untouched.
    - Output: reprocess_<original_name>.dat

    Returns output path as a string.
    """
    import re
    from pathlib import Path

    input_path = Path(input_path)
    if input_path.suffix.lower() != ".dat":
        raise ValueError("Expected a .dat input file")

    out_path = input_path.with_name("reprocess_" + input_path.name)

    stoich_re = re.compile(r"^(\d+)([A-Z][A-Za-z0-9()]*)$")

    def expand_side(side: str) -> str:
        terms = [t.strip() for t in side.split("+")]
        expanded = []

        for term in terms:
            if not term:
                continue

            m = stoich_re.match(term)
            if m:
                n = int(m.group(1))
                sp = m.group(2)
                expanded.extend([sp] * n)
            else:
                expanded.append(term)

        return " + ".join(expanded)

    with open(input_path, "r") as fin, open(out_path, "w") as fout:
        for line in fin:
            chunk = line[6:70].rstrip("\n")
            chunk = chunk.strip()
            if not chunk:
                continue

            # normalize separators
            chunk = re.sub(r"\s*\+\s*", " + ", chunk)
            chunk = re.sub(r"\s*=\s*", " = ", chunk)

            if " = " not in chunk:
                continue

            lhs, rhs = chunk.split(" = ", 1)

            lhs = expand_side(lhs)
            rhs = expand_side(rhs)

            fout.write(f"{lhs} = {rhs}\n")
    print(f"Done reprocessing {input_path}\n")
    return str(out_path)


def reprocess_stand_file(input_path):
    """

    - Skip the FIRST line of the file.
    - For all other lines: take only the first 50 columns, then apply:

    1) Remove ALL underscores everywhere.
    2) Replace "Oxyrane" (case-insensitive) with C2H4O.
    3) Replace '=>' with '='.
    4) Keep parentheses (never remove '(' or ')').
    5) Remove excitation carets inside/around parentheses:
         O_2(a^1Delta__g^) -> O2(a1Deltag)
         A(X^Y) -> A(XY)
         A(anything^) -> A(anything)
       (i.e., remove '^' characters; parentheses remain)
    6) Star states:
         A^^*^, A^*^, A^* (token-bounded) -> A(*)
         e.g., CH2^* -> CH2(*)
    7) Charge cleanup:
         - If + or - is part of a token (not a separator), canonicalize to ^+ or ^-
           Examples:
             A^+^^ -> A^+
             Cl+   -> Cl^+
             O2-   -> O2^-
           This does NOT touch separator '+' with spaces around it.
         - Ensures no space between species and its ^+/^-
    8) Replace 'gamma' (case sensitive) with 'HV'.
    9) Normalize electron notations to 'e^-':
         e-, e^-, ^e-^, e^-^^, etc. -> e^-
         Only when 'e' is standalone (bounded by start/space/^ on left).
    10) Replace exact 'h{nu}' with 'HV'.

    Output: reprocess_<original_name>.dat
    Returns output path as a string.
    """
    import re
    from pathlib import Path

    input_path = Path(input_path)
    if input_path.suffix.lower() != ".dat":
        raise ValueError("Expected a .dat input file")

    out_path = input_path.with_name("reprocess_" + input_path.name)

    # Electron patterns (standalone e with various notations)
    e_re = re.compile(r"(?:(?<=\s)|(?<=\^)|^)\^?e(?:\^\-|\-|\^\-?\^*|\-+\^*)", re.IGNORECASE)

    # Star states inside a token: CH2^*, CH2^^*^, etc.
    star_re = re.compile(r"^([A-Za-z0-9][A-Za-z0-9()]*)\^{0,2}\*\^{0,}$")

    # Collapse caret-junk around charge: ^^+^^ -> ^+  (and same for -)
    caret_charge_junk_re = re.compile(r"\^{1,}([+-])\^{0,}")

    # Plain token charges: Cl+ / O2- -> Cl^+ / O2^- (must be token-wise)
    plain_charge_tok_re = re.compile(r"^([A-Za-z0-9][A-Za-z0-9()]*)\s*([+-])$")

    # Remove any caret that is not part of charge (we'll protect ^+ ^- and e^- first)
    caret_noncharge_re = re.compile(r"\^(?![+-])")

    def normalize_tokens(s: str) -> str:
        toks = s.split()  # separator '+' becomes its own token if spaced; good
        out = []

        for t in toks:
            # Normalize electron-like tokens early if they appear as a token
            if t.lower().startswith("e"):
                # token-wise normalization is handled by regex on whole string later
                out.append(t)
                continue

            # Star state token
            m_star = star_re.match(t)
            if m_star:
                out.append(m_star.group(1) + "(*)")
                continue

            # Collapse caret-junk charge inside token (e.g. A^+^^)
            t2 = caret_charge_junk_re.sub(r"^\1", t)

            # Plain charge token (e.g. Cl+)
            m_plain = plain_charge_tok_re.match(t2)
            if m_plain:
                base, sign = m_plain.group(1), m_plain.group(2)
                out.append(f"{base}^{sign}")
            else:
                out.append(t2)

        return " ".join(out)

    with open(input_path, "r") as fin, open(out_path, "w") as fout:
        _ = next(fin, None)  # skip first line

        for line in fin:
            s = line[:50].rstrip("\n")
            if not s.strip():
                continue

            # (10) h{nu} -> HV (exact)
            s = s.replace("h{nu}", "HV")

            # (2) Oxyrane -> C2H4O (case-insensitive)
            s = re.sub(r"(?i)\boxyrane\b", "C2H4O", s)

            # (3) => -> =
            s = s.replace("=>", "=")

            # (8) gamma -> HV (case sensitive)
            s = s.replace("gamma", "HV")

            # (1) remove underscores everywhere
            s = s.replace("_", "")

            # (9) normalize electrons -> e^-
            s = e_re.sub("e^-", s)

            # Token-wise star + charge processing (rules 6 & 7)
            s = normalize_tokens(s)

            # Protect electron so its caret isn't removed below
            s = s.replace("e^-", "__ELECTRON__")

            # (5) remove remaining carets used for excited-state notation, keep parentheses
            # This turns O(^1^D) -> O(1D), O2(a^1Deltag^) -> O2(a1Deltag)
            s = caret_noncharge_re.sub("", s)

            # Restore electron
            s = s.replace("__ELECTRON__", "e^-")

            # Final cleanup spacing around separators
            s = re.sub(r"\s*=\s*", " = ", s)
            s = re.sub(r"\s+", " ", s).strip()

            fout.write(s + "\n")
    print(f"Done reprocessing {input_path}\n")
    return str(out_path)


