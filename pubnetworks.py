import re
import os
from collections import Counter, OrderedDict
from pathlib import Path
from inputchem_backend import write_unique_species_list, sort_rxns_reactant

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

# legacy
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

# https://doi.org/10.1038/s41586-023-05902-2 , data presented here
def moses_vII(path, allowed_elements, metadata="N", info = "N"):
    """
    Master function for Moses chemical network.
    path: folder containing ReadMe and input_MosesII .txt files
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
        #line = re.sub(r'^\s*\d+\)\s*', '', line)

        # cut off kinetics (everything starting with k, k0, koo)
        line = line[:60].rstrip()

        if '=' not in line:
            return[]

        lhs, rhs = line.split('=', 1)

        parts = lhs.split('+') + rhs.split('+')
        return [p.strip() for p in parts if p.strip()]


    # ---------------- File Setup ----------------
    tag = "_".join(sorted(el for el in allowed_elements))
    parent_folder = os.path.dirname(path)
    output_file = os.path.join(parent_folder, 
                               f"MOSESII_output_{tag}.dat")
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
    txt_file = os.path.join(path, "input_MosesII.txt")
    #print(txt_file)
    if not os.path.exists(txt_file):
        txt_file = os.path.join(path, txt_files[0])  # fallback
        print(f"Specified file not found, using {txt_file}")

    # ---------------- Process File ----------------
    stop_line_pattern = " STOP_EOF"
    reactions_started = False

    with open(txt_file, "r") as f:
        for lineno, line in enumerate(f, 1):
            raw = line.rstrip("\n")
            
            # First 70 characters as header
            if info == 'Y':
                if lineno == 1:
                    kept.append("\n\n\n")
                    kept.append(raw[:70])
                    continue


            # Stop condition
            if stop_line_pattern in raw or lineno > 1450:
                break

            # Skip until "REACTIONS:" is found
            if not reactions_started:
                if info == 'Y':
                    kept.append(raw)
                if "STOP SPECIES" in raw:
                    reactions_started = True
                continue

            # Extract reaction line
            if re.match(r'^\s*[A-Z0-9]', raw):
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
    out_file = os.path.join(os.path.dirname(path), f"velliet_output_{tag}.dat")

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


THIS PART IS BROKEN DOWN INTO TWO SUB-PARTS:

A) CHEMICAL REACTION FILE REPROCESSING: To a homogenous nomeclature of the 
reaction networks. Standardization of all different reaction writing proc-
dures used by different authors.

INPUT: The CHEMICAL NETWORK FILES produced in PART I (without metadata!!!).
OUTPUT: Files with name "reprocess_{input_file_name}.dat".
These files contain the reactions in a one standardized uniform format.

IMP :: GENERAL CHARACTERISTICS OF THE FILES PRODUCED
- HV = Used to denote light (h\nu) across all the files.
- BRACKETS () = Used to write the excited states. 
                e.g: O(3P), O(1S), O(X2Pig), CH2(*) etc
- ^+ or ^- or ^-- = In general the caret ^ is used to denote ions in the file.

- FOR STAND network file processing only; OXYRANE (in input file) = C2H4O (its
chemical formula) in the output file

- For VELLIET_VENOT; O3P(in input file) is written as O(in output file). 
- FOR VV; " OOH " changes to HO2.
- For Hu; we add " + T" to lhs for the thermal decomposition reactions

"""


def reprocess_velliet_file(path, rate_const="Y"):
    """
    Reprocess a Velliet/Venot file into "A + B = C + D" reactions.

    Existing behavior (unchanged):
    - Skip blank lines
    - Skip lines starting with lowercase OR starting with '-'
    - Read only first 110 columns
    - Convert space-blocks:
        3..11 spaces  -> " + "
        >=30 spaces   -> " = "
        otherwise     -> " "
    - Cleanup doubles and normalize whitespace
    - Replace token O3P -> O

    NEW behavior (your request):
    - Some sections contain 3-body reactions and must get " + M" on BOTH sides.
      Triggered by section header lines in the file:
        Start adding M at: reactions_1.dat
          Stop adding M at: reactions_3.dat if present, else at reactions_4.dat
        Start adding M at: reactions_5.dat
          Stop adding M at: reactions_7.dat
        Start adding M at: reactions_9.dat
          Stop adding M at: reactions_10.dat

      Applies only to reaction lines (same filtering rule: not starting with lowercase, not starting with '-').

    RATE CONSTANT EXTRA OUTPUT (your request):
    - If rate_const == "Y":
        also write a second file that keeps the processed reaction, then appends
        the original line content AFTER column 110 (the part we previously ignored).
      Output name:
        reprocess_rate_const_<original_filename>
    """
    path = Path(path)
    out_path = path.with_name("reprocess_" + path.name)

    rate_const = str(rate_const).strip().upper()
    out_path_rc = path.with_name("reprocess_rate_const_" + path.name) if rate_const == "Y" else None

    processed = []
    processed_rc = []  # only used if rate_const == "Y"

    # --- section control for "+ M" injection ---
    add_M = False
    stop_after_1_2 = None  # will be decided on-the-fly: 3 if encountered, else stop at 4

    # robust header detector
    hdr_re = re.compile(r"\breactions_(\d+)\.dat\b", re.IGNORECASE)

    def add_M_both_sides(reaction: str) -> str:
        # add " + M" to both sides of " = " exactly once
        if " = " not in reaction:
            return reaction

        lhs, rhs = reaction.split(" = ", 1)

        # If M already present as a token, do nothing for that side
        lhs_tokens = lhs.split(" + ")
        rhs_tokens = rhs.split(" + ")

        if "M" not in lhs_tokens:
            lhs = lhs.rstrip() + " + M"
        if "M" not in rhs_tokens:
            rhs = rhs.rstrip() + " + M"

        return lhs + " = " + rhs

    with open(path) as f:
        for line in f:
            if not line.strip():
                continue

            # keep full original (no newline) for rate-constant output
            original_full = line.rstrip("\n")

            # ---- NEW: detect section headers BEFORE skipping lowercase lines ----
            m_hdr = hdr_re.search(line)
            if m_hdr:
                n = int(m_hdr.group(1))

                if n == 1:
                    add_M = True
                    stop_after_1_2 = None
                elif n == 2:
                    add_M = True
                elif n == 3:
                    add_M = False
                    stop_after_1_2 = 3
                elif n == 4:
                    # stop after reactions_1/reactions_2 if reactions_3 never appeared
                    if stop_after_1_2 is None:
                        add_M = False
                        stop_after_1_2 = 4
                elif n == 5:
                    add_M = True
                elif n == 7:
                    add_M = False
                elif n == 9:
                    add_M = True
                elif 10 <= n <= 14:
                    add_M = False

                # header lines are not reactions
                continue

            # ---- existing filtering (unchanged) ----
            if line[0].islower() or line.startswith("-"):
                continue

            # the part you already used for reaction parsing
            line = line[:110].rstrip()

            # split into text and space blocks
            parts = re.split(r"(\s+)", line)

            new = []
            for p in parts:
                if p.isspace():
                    n = len(p)
                    if n >= 3 and n <= 11:
                        new.append(" + ")
                    elif n >= 30:
                        new.append(" = ")
                    else:
                        new.append(" ")
                else:
                    new.append(p.strip())

            reaction = "".join(new)

            # replace O3P token with O
            reaction = re.sub(r"\bO3P\b", "O", reaction)
            reaction = re.sub(r"\bOOH\b", "HO2", reaction)

            # clean up accidental doubles
            reaction = re.sub(r"\s+\+\s+\+", " + ", reaction)
            reaction = re.sub(r"\s+=\s+=", " = ", reaction)
            reaction = re.sub(r"\s+", " ", reaction).strip()

            # ---- NEW: add "+ M" for the specified sections ----
            if add_M:
                reaction = add_M_both_sides(reaction)

                # keep formatting clean after adding M
                reaction = re.sub(r"\s+", " ", reaction).strip()
                reaction = re.sub(r"\s*\+\s*", " + ", reaction)
                reaction = re.sub(r"\s*=\s*", " = ", reaction)

            processed.append(reaction)

            # ---- NEW: optional rate-constant output (append ignored tail after col 110) ----
            if rate_const == "Y":
                tail = ""
                if len(original_full) > 110:
                    tail = original_full[110:]  # keep EXACTLY what was previously ignored
                tail = tail.rstrip()

                if tail:
                    processed_rc.append(reaction.ljust(50)[:50] + "\t" + tail)
                else:
                    processed_rc.append(reaction.ljust(50)[:50])

    with open(out_path, "w") as f:
        for r in processed:
            f.write(r + "\n")

    if rate_const == "Y":
        with open(out_path_rc, "w") as f:
            for r in processed_rc:
                f.write(r + "\n")

    print(f"Done reprocessing {path}\n")
    if rate_const == "Y":
        print(f"Also wrote rate-constant version: {out_path_rc}\n")

    return out_path


def reprocess_agundez_file(input_path, rate_const="Y"):
    """
    Reprocess an Agundez_output_*.dat file:
      - keep only first 50 characters of each line
      - convert tokens like A_BC -> A(BC), where BC goes until next whitespace
      - output is exactly 50 characters wide per line (padded or trimmed)
      - writes to reprocess_<original_name>.dat

    RATE CONSTANT EXTRA OUTPUT (your request):
    - If rate_const == "Y":
        also write a second file that keeps the processed 50-col chunk, then appends
        the original line content AFTER column 50 (the part we previously ignored).
      Output name:
        reprocess_rate_const_<original_filename>.dat
    """
    input_path = Path(input_path)
    if input_path.suffix.lower() != ".dat":
        raise ValueError("Expected a .dat input file")

    out_path = input_path.with_name("reprocess_" + input_path.name)

    rate_const = str(rate_const).strip().upper()
    out_path_rc = (
        input_path.with_name("reprocess_rate_const_" + input_path.name)
        if rate_const == "Y"
        else None
    )

    # Match a "word" chunk containing underscore, stopping at whitespace.
    # Example matches: A_BC, A_3X, O_1D, CH3_foo
    pattern = re.compile(r"([^\s_]+)_([^\s]+)")

    with open(input_path, "r") as fin, open(out_path, "w") as fout:
        fout_rc = open(out_path_rc, "w") if rate_const == "Y" else None
        try:
            for line in fin:
                original_full = line.rstrip("\n")

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

                # Optional: append ignored tail after col 50
                if fout_rc is not None:
                    tail = ""
                    if len(original_full) > 50:
                        tail = original_full[50:]
                    tail = tail.rstrip()

                    if tail:
                        fout_rc.write(chunk.ljust(50)[:50] + "\t" + tail + "\n")
                    else:
                        fout_rc.write(chunk.ljust(50)[:50] + "\n")
        finally:
            if fout_rc is not None:
                fout_rc.close()

    print(f"Done reprocessing {input_path}\n")
    if rate_const == "Y":
        print(f"Also wrote rate-constant version: {out_path_rc}\n")
    return str(out_path)

def reprocess_hu_file(input_path, rate_const="Y"):

    input_path = Path(input_path)
    if input_path.suffix.lower() != ".dat":
        raise ValueError("Expected a .dat input file")

    out_path = input_path.with_name("reprocess_" + input_path.name)

    rate_const = str(rate_const).strip().upper()
    out_path_rc = (
        input_path.with_name("reprocess_rate_const_" + input_path.name)
        if rate_const == "Y"
        else None
    )

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
        fout_rc = open(out_path_rc, "w") if rate_const == "Y" else None
        try:
            for line in fin:
                original_full = line.rstrip("\n")

                if not line.startswith(("R", "M", "T")):
                    continue

                starts_with_M = line.startswith("M")
                starts_with_T = line.startswith("T")

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
                if starts_with_T:
                    lhs = lhs.strip()
                    if lhs:
                        lhs += " + T"

                out_line = f"{lhs} = {rhs}".strip()
                out_line = normalize_caret_notation_in_line(out_line)

                fout.write(out_line + "\n")

                # Optional: append ignored tail after the chunk slice end (index 37)
                if fout_rc is not None:
                    tail = ""
                    if len(original_full) > 37:
                        tail = original_full[37:]  # EXACTLY what was ignored
                    tail = tail.rstrip()

                    if tail:
                        fout_rc.write(out_line.ljust(50)[:50] + "\t" + tail + "\n")
                    else:
                        fout_rc.write(out_line.ljust(50)[:50]  + "\n")
        finally:
            if fout_rc is not None:
                fout_rc.close()

    print(f"Done reprocessing {input_path}\n")
    if rate_const == "Y":
        print(f"Also wrote rate-constant version: {out_path_rc}\n")
    return str(out_path)

def reprocess_vulcan_file(input_path, rate_const="Y"):
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

    RATE CONSTANT EXTRA OUTPUT (your request):
    - If rate_const == "Y":
        also write a second file that keeps the processed out_line, then appends
        the original line content AFTER column 42 (the part we previously ignored).
      Output name:
        reprocess_rate_const_<original_filename>

    Returns output path as a string.
    """
    input_path = Path(input_path)
    if input_path.suffix.lower() != ".dat":
        raise ValueError("Expected a .dat input file")

    out_path = input_path.with_name("reprocess_" + input_path.name)

    rate_const = str(rate_const).strip().upper()
    out_path_rc = (
        input_path.with_name("reprocess_rate_const_" + input_path.name)
        if rate_const == "Y"
        else None
    )

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
        fout_rc = open(out_path_rc, "w") if rate_const == "Y" else None
        try:
            for line in fin:
                original_full = line.rstrip("\n")

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

                # Optional: append ignored tail after slice end index 42
                if fout_rc is not None:
                    tail = ""
                    if len(original_full) > 42:
                        tail = original_full[44:]  # EXACTLY what was ignored
                    tail = tail.rstrip()

                    if tail:
                        fout_rc.write(out_line.ljust(50)[:50] + "\t" + tail + "\n")
                    else:
                        fout_rc.write(out_line.ljust(50)[:50] + "\n")
        finally:
            if fout_rc is not None:
                fout_rc.close()

    print(f"Done reprocessing {input_path}\n")
    if rate_const == "Y":
        print(f"Also wrote rate-constant version: {out_path_rc}\n")
    return str(out_path)

# Legacy code for original MOSES 2011
def reprocess_moses_file(input_path, rate_const="Y"):
    """
    Moses reprocess:

    - Read columns 7–70 (1-based) => Python slice [6:70]
    - Normalize spacing around '+' and '='
    - Expand stoichiometric prefixes like 2H, 2OH, 2HO2, 3O, etc.
      ONLY when the number is a PREFIX and the species starts with a capital letter.
      Molecules like H2O are untouched.
    - Output: reprocess_<original_name>.dat

    RATE CONSTANT EXTRA OUTPUT (your request):
    - If rate_const == "Y":
        also write a second file that keeps the processed reaction, then appends
        the original line content AFTER column 70 (the part we previously ignored).
      Output name:
        reprocess_rate_const_<original_filename>

    Returns output path as a string.
    """

    input_path = Path(input_path)
    if input_path.suffix.lower() != ".dat":
        raise ValueError("Expected a .dat input file")

    out_path = input_path.with_name("reprocess_" + input_path.name)

    rate_const = str(rate_const).strip().upper()
    out_path_rc = (
        input_path.with_name("reprocess_rate_const_" + input_path.name)
        if rate_const == "Y"
        else None
    )

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
        fout_rc = open(out_path_rc, "w") if rate_const == "Y" else None
        try:
            for line in fin:
                original_full = line.rstrip("\n")

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

                out_line = f"{lhs} = {rhs}"
                fout.write(out_line + "\n")

                # Optional: append ignored tail after slice end index 70
                if fout_rc is not None:
                    tail = ""
                    if len(original_full) > 70:
                        tail = original_full[70:]  # EXACTLY what was ignored
                    tail = tail.rstrip()

                    if tail:
                        fout_rc.write(out_line.ljust(50)[:50] + "\t" + tail + "\n")
                    else:
                        fout_rc.write(out_line.ljust(50)[:50]  + "\n")
        finally:
            if fout_rc is not None:
                fout_rc.close()

    print(f"Done reprocessing {input_path}\n")
    if rate_const == "Y":
        print(f"Also wrote rate-constant version: {out_path_rc}\n")
    return str(out_path)

# code for MOSES network as presented in
# https://doi.org/10.1038/s41586-023-05902-2 , by SHAMI
def reprocess_mosesvii_file(input_path, rate_const="Y"):
    """
    Moses reprocess:

    - Read columns 0–60 (1-based) => Python slice [0:59]
    - Normalize spacing around '+' and '='
    - Expand stoichiometric prefixes like 2H, 2OH, 2HO2, 3O, etc.
      ONLY when the number is a PREFIX and the species starts with a capital letter.
      Molecules like H2O are untouched.
    - Output: reprocess_<original_name>.dat

    RATE CONSTANT EXTRA OUTPUT (your request):
    - If rate_const == "Y":
        also write a second file that keeps the processed reaction, then appends
        the original line content AFTER column 60 (the part we previously ignored).
      Output name:
        reprocess_rate_const_<original_filename>

    Returns output path as a string.
    """

    input_path = Path(input_path)
    if input_path.suffix.lower() != ".dat":
        raise ValueError("Expected a .dat input file")

    out_path = input_path.with_name("reprocess_" + input_path.name)

    rate_const = str(rate_const).strip().upper()
    out_path_rc = (
        input_path.with_name("reprocess_rate_const_" + input_path.name)
        if rate_const == "Y"
        else None
    )

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
        fout_rc = open(out_path_rc, "w") if rate_const == "Y" else None
        try:
            for line in fin:
                original_full = line.rstrip("\n")

                chunk = line[0:59].rstrip("\n")
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

                out_line = f"{lhs} = {rhs}"
                fout.write(out_line + "\n")

                # Optional: append ignored tail after slice end index 60
                if fout_rc is not None:
                    tail = ""
                    if len(original_full) > 59:
                        tail = original_full[59:]  # EXACTLY what was ignored
                    tail = tail.rstrip()

                    if tail:
                        fout_rc.write(out_line.ljust(50)[:50] + "\t" + tail + "\n")
                    else:
                        fout_rc.write(out_line.ljust(50)[:50]  + "\n")
        finally:
            if fout_rc is not None:
                fout_rc.close()

    print(f"Done reprocessing {input_path}\n")
    if rate_const == "Y":
        print(f"Also wrote rate-constant version: {out_path_rc}\n")
    return str(out_path)



def reprocess_stand_file(input_path, rate_const="Y"):
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

    RATE CONSTANT EXTRA OUTPUT (your request):
    - If rate_const == "Y":
        also write a second file that keeps the processed line, then appends
        the original line content AFTER column 50 (the part we previously ignored).
      Output name:
        reprocess_rate_const_<original_filename>
    """
    input_path = Path(input_path)
    if input_path.suffix.lower() != ".dat":
        raise ValueError("Expected a .dat input file")

    out_path = input_path.with_name("reprocess_" + input_path.name)

    rate_const = str(rate_const).strip().upper()
    out_path_rc = (
        input_path.with_name("reprocess_rate_const_" + input_path.name)
        if rate_const == "Y"
        else None
    )

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
        fout_rc = open(out_path_rc, "w") if rate_const == "Y" else None
        try:
            _ = next(fin, None)  # skip first line

            for line in fin:
                original_full = line.rstrip("\n")

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
                s = s.replace("a1Deltag", "1Delta")
                s = re.sub(r"\s*=\s*", " = ", s)
                s = re.sub(r"\s+", " ", s).strip()

                fout.write(s + "\n")

                # Optional: append ignored tail after col 50
                if fout_rc is not None:
                    tail = ""
                    if len(original_full) > 50:
                        tail = original_full[50:]  # EXACTLY what was ignored
                    tail = tail.rstrip()

                    if tail:
                        fout_rc.write(s.ljust(50)[:50] + "\t" + tail + "\n")
                    else:
                        fout_rc.write(s.ljust(50)[:50] + "\n")
        finally:
            if fout_rc is not None:
                fout_rc.close()

    print(f"Done reprocessing {input_path}\n")
    if rate_const == "Y":
        print(f"Also wrote rate-constant version: {out_path_rc}\n")
    return str(out_path)


########## Antonio 2007 chemistry.dat processing       #########



#================================= UNIQUE REACTION PRODUCTION =========================================

def build_unique_reactions(reprocess_paths, long_table="N"):
    """
    Build unique-reaction tables across 6(or n [can be expanded]) *reprocessed* network files.

    Inputs
    ------
    reprocess_paths : list[str] | tuple[str]
        Paths to the 6 reprocessed .dat files (Hu, Agundez, Velliet, Moses, VULCAN, Stand),
        in ANY order. The function derives each network name from the filename:
            reprocess_<NETWORK>_output_<ELEMENTS>.dat
        and derives the <ELEMENTS> tag from between 'output_' and '.dat'.
    long_table : str
        "Y" to also write the long-table file, otherwise "N".

    Outputs
    -------
    1) Uniq_reactions_<ELEMENTS>.dat
        One reaction per row (unique across all inputs).
    2) Uniq_reactions_metadata_<ELEMENTS>.dat
        First column (70 chars) is the reaction text.
        Then 6 columns (10 chars each) giving the line number (0001 format) in each network,
        or 0000 if absent.
    3) (optional) Uniq_reactions_longtable_<ELEMENTS>.dat
        6 columns (one per network). Each row shows the reaction as written in that network.
        Blank if absent.

    Canonicalization rules for matching
    ----------------------------------
    - Reaction separator is ' = ' (space-equals-space).
    - Species separator is ' + ' (space-plus-space).
      This avoids splitting ion charges like '^+' or '++' since those are not space-bounded.
    - Ignore standalone photon token: exactly 'HV' as a token (i.e., bounded by spaces via split).
    - Ignore only the literal parentheses characters '(' and ')', but keep their contents:
        O(1D) -> O1D   (for matching only)
    - Order within each side does not matter:
        A + BC = AB + C matches BC + A = C + AB
    - Direction is preserved (lhs != rhs). If you want to treat reverse as identical, say so.

    Notes
    -----
    - Assumes each input line is already "reprocessed" and reaction-like.
    - The first occurrence per file is used for metadata (line number) and long-table text.
    """

    nof = len(reprocess_paths)
    if not isinstance(reprocess_paths, (list, tuple)) or len(reprocess_paths) != 6:
        print(f"\n {nof} files given \n reprocess_paths must be a list/tuple of exactly 6 file paths.")
        print("Handles more than 6 files as well. Just a sanity check....")

    long_table = (long_table or "N").strip().upper()
    if long_table not in {"Y", "N"}:
        long_table = "N"

    # -------------------------
    # helpers: name + tag
    # -------------------------
    def network_name_from_path(p):
        base = os.path.basename(p)
        m = re.search(r"^reprocess_(.+?)_output_", base)
        if not m:
            print("Something wrong with the file names.\n Check if the network name is present in filename")
            # fallback: try without reprocess_ prefix
            m = re.search(r"^(.+?)_output_", base)
        return m.group(1) if m else os.path.splitext(base)[0]

    def tag_from_any_path(p):
        base = os.path.basename(p)
        m = re.search(r"output_(.+?)\.dat$", base)
        return m.group(1) if m else "UNKNOWN"

    tag = None
    for p in reprocess_paths:
        t = tag_from_any_path(p)
        if t != "UNKNOWN":
            tag = t
            break
    if tag is None:
        tag = "UNKNOWN"

    # -------------------------
    # parsing / canonicalization
    # -------------------------
    SEP_RXN = " = "
    SEP_SPECIES = " + "


    def strip_paren_chars(species):
    # keep contents, drop only parentheses characters
        return species.replace("(", "").replace(")", "")
 
    def normalize_species_for_matching(sp):
        sp2 = strip_paren_chars(sp)

        # NEW: normalize OH/HO with charges, assuming at most one caret in the token
        # Accept: OH^+, OH^++, HO^---, OH++, HO-
        # Preserve multiplicity: OH^+ != OH^++
        # Canonicalize base: OH/HO -> HO
        m = re.fullmatch(r"([A-Z]{2})(?:\^)?([+-]+)", sp2)
        if m:
            base, charge = m.group(1), m.group(2)
            if base in {"OH", "HO"}:
                # reject mixed + and - in the same token
                if "+" in charge and "-" in charge:
                    return sp2
                return "".join(sorted(base)) + "^" + charge  # e.g., HO^++, HO^-

        # Existing rule: reorder only if token is exactly [A-Z]+
        if re.fullmatch(r"[A-Z]+", sp2):
            sp2 = "".join(sorted(sp2))

        return sp2



    def canonicalize_side(side_text):
        # split by ' + ' exactly; do NOT split on '+' used as charge
        parts = [p.strip() for p in side_text.split(SEP_SPECIES)]
        out = []
        for sp in parts:
            if not sp:
                continue
            # ignore photon only if it's the standalone token HV
            if sp == "HV":
                continue
            sp2 = normalize_species_for_matching(sp)
            out.append(sp2)
        return tuple(sorted(out))

    def canonical_key(reaction_line):
        # returns (lhs_tuple_sorted, rhs_tuple_sorted) or None
        if SEP_RXN not in reaction_line:
            return None
        lhs, rhs = reaction_line.split(SEP_RXN, 1)
        lhs_t = canonicalize_side(lhs.strip())
        rhs_t = canonicalize_side(rhs.strip())
        # If everything was stripped (e.g., only HV), ignore
        if not lhs_t and not rhs_t:
            return None
        return (lhs_t, rhs_t)

    # -------------------------
    # ingest
    # -------------------------
    net_order = [network_name_from_path(p) for p in reprocess_paths]
    print(net_order)

    # stable order (to keep outputs deterministic):
    # If the "usual" names exist, use that order, otherwise use discovered order.
    preferred = ["Hu", "Agundez", "VULCAN", "velliet", "MOSESII", "stand"]

    def rank(n):
        return preferred.index(n) if n in preferred else 10_000 + net_order.index(n)
    net_order = sorted(net_order, key=rank)

    # map network name -> path
    net_path = {network_name_from_path(p): p for p in reprocess_paths}

    # database:
    # key -> {
    #   "repr": first_seen_reaction_string,
    #   "line": {net: line_number_int or 0},
    #   "text": {net: reaction_string or ""}
    # }
    db = OrderedDict()

    for net in net_order:
        path = net_path[net]
        with open(path, "r") as f:
            for ln, raw in enumerate(f, 1):
                line = raw.rstrip("\n").strip()
                if not line:
                    continue
                # enforce consistent separator expectation
                if SEP_RXN not in line:
                    continue

                key = canonical_key(line)
                if key is None:
                    continue

                if key not in db:
                    db[key] = {
                        "repr": line,  # keep original formatting from first occurrence globally
                        "line": {n: 0 for n in net_order},
                        "text": {n: "" for n in net_order},
                    }

                # first occurrence per network
                if db[key]["line"][net] == 0:
                    db[key]["line"][net] = ln
                    db[key]["text"][net] = line

    # -------------------------
    # write outputs
    # -------------------------
    out1 = f"uniq_reactions_{tag}.dat"
    out2 = f"uniq_reactions_metadata_{tag}.dat"
    out3 = f"uniq_reactions_longtable_{tag}.dat"

    # --- file 1: unique reactions list
    with open(out1, "w") as f1:
        for entry in db.values():
            f1.write(entry["repr"] + "\n")

    # --- file 2: metadata (fixed widths)
    RXN_W = 50
    COL_W = 10

    def fmt_rxn(s):
        s = (s or "").strip()
        return (s[:RXN_W]).ljust(RXN_W)

    def fmt_ln(n):
        try:
            n = int(n)
        except Exception:
            n = 0
        if n < 0:
            n = 0
        if n > 9999:
            n = 9999
        return f"{n:04d}".ljust(COL_W)

    header = [("Reaction".ljust(RXN_W))] + [net.ljust(COL_W) for net in net_order]

    with open(out2, "w") as f2:
        f2.write("\t".join(header) + "\n")
        for entry in db.values():
            row = [fmt_rxn(entry["repr"])]
            row.extend(fmt_ln(entry["line"][net]) for net in net_order)
            f2.write("\t".join(row) + "\n")

   
    # --- file 3: optional long table (FIXED 50-COLUMN WIDTH + COLUMN TOTALS)
    if long_table == "Y":
        COL_W_LONG = 50

        def fmt_long(s):
            s = (s or "").strip()
            return s[:COL_W_LONG].ljust(COL_W_LONG)

        # initialize counters
        counts = {net: 0 for net in net_order}

        with open(out3, "w") as f3:
            # header
            f3.write("".join(fmt_long(net) for net in net_order) + "\n")

            # rows
            for entry in db.values():
                row_parts = []
                for net in net_order:
                    txt = entry["text"][net]
                    if txt:
                        counts[net] += 1
                    row_parts.append(fmt_long(txt))
                f3.write("".join(row_parts) + "\n")

            # final totals row
            f3.write("".join(fmt_long(str(counts[net])) for net in net_order) + "\n")

    print("done finding unique reactions....")
    return {
        "tag": tag,
        "networks": net_order,
        "out_unique": out1,
        "out_metadata": out2,
        "out_longtable": out3 if long_table == "Y" else None,
        "n_unique": len(db),
    }



#================================= UNIQUE REACTION PRODUCTION - VII =========================================
# SAME AS ABOVE JUST HANDLES THE RATE CONSTANT AS WELL
# Gives the whole line as output after unique reaction.

def build_unique_rxns_ratek(reprocess_rateconst_paths, long_table="N"):
    """
    Like build_unique_reactions(), but for reprocess_rate_constant*.dat files where:
      - reaction text is in the first 50 chars (fixed-width)
      - tail text is line[50:] (original ignored chunk)

    Outputs:
    1) uniq_reactions_rateconst_<TAG>.dat
       One unique reaction per row (uses first-seen reaction text, stripped).
    2) uniq_reactions_rateconst_metadata_<TAG>.dat
       First line: Reaction (50 chars) + 7 columns of line numbers (0001 format) per network.
       Then, UNDER EACH reaction row, prints the tails in this exact order:
         Hu, Agundez, VULCAN, velliet, MOSES, stand, agm2007
       Format:
         NET \\t :: \\t TAIL
       If missing or tail empty: NET \\t :: \\t !!!!!
    3) (optional) uniq_reactions_rateconst_longtable_<TAG>.dat
       Same idea as your long table: one column per network showing reaction text (first 50 chars).
    """

    nof = len(reprocess_rateconst_paths)
    print(reprocess_rateconst_paths)
    if not isinstance(reprocess_rateconst_paths, (list, tuple)) or len(reprocess_rateconst_paths) < 1:
        raise ValueError("reprocess_rateconst_paths must be a list/tuple of file paths.")

    long_table = (long_table or "N").strip().upper()
    if long_table not in {"Y", "N"}:
        long_table = "N"

    # -------------------------
    # helpers: name + tag
    # -------------------------
    def network_name_from_path(p):
        base = os.path.basename(p)

        # Try both:
        # reprocess_rate_constant<NETWORK>_output_...
        # reprocess_rate_constant<something> (fallbacks)
        m = re.search(r"^reprocess_rate_const_(.+?)_output_", base)
        if not m:
            # If your file is actually named reprocess_rate_constant<original> with no _output_ pattern,
            # fall back to your old extractor.
            m = re.search(r"^reprocess_(.+?)_output_", base)
            if not m:
                m = re.search(r"^(.+?)_output_", base)

        return m.group(1) if m else os.path.splitext(base)[0]

    def tag_from_any_path(p):
        base = os.path.basename(p)
        m = re.search(r"output_(.+?)\.dat$", base)
        return m.group(1) if m else "UNKNOWN"

    tag = None
    for p in reprocess_rateconst_paths:
        t = tag_from_any_path(p)
        if t != "UNKNOWN":
            tag = t
            break
    if tag is None:
        tag = "UNKNOWN"

    # -------------------------
    # parsing / canonicalization (MATCHING USES ONLY FIRST 50 CHARS)
    # -------------------------
    SEP_RXN = " = "
    SEP_SPECIES = " + "

    def strip_paren_chars(species):
        return species.replace("(", "").replace(")", "")

    def normalize_species_for_matching(sp):
        sp2 = strip_paren_chars(sp)

        # OH/HO charge normalization (your existing behavior)
        m = re.fullmatch(r"([A-Z]{2})(?:\^)?([+-]+)", sp2)
        if m:
            base, charge = m.group(1), m.group(2)
            if base in {"OH", "HO"}:
                if "+" in charge and "-" in charge:
                    return sp2
                return "".join(sorted(base)) + "^" + charge

        if re.fullmatch(r"[A-Z]+", sp2):
            sp2 = "".join(sorted(sp2))

        return sp2

    def canonicalize_side(side_text):
        parts = [p.strip() for p in side_text.split(SEP_SPECIES)]
        out = []
        for sp in parts:
            if not sp:
                continue
            if sp == "HV":
                continue
            out.append(normalize_species_for_matching(sp))
        return tuple(sorted(out))

    def canonical_key(reaction_line):
        if SEP_RXN not in reaction_line:
            return None
        lhs, rhs = reaction_line.split(SEP_RXN, 1)
        lhs_t = canonicalize_side(lhs.strip())
        rhs_t = canonicalize_side(rhs.strip())
        if not lhs_t and not rhs_t:
            return None
        return (lhs_t, rhs_t)

    # -------------------------
    # ingest
    # -------------------------
    discovered = [network_name_from_path(p) for p in reprocess_rateconst_paths]

    # Your requested order for tail blocks
    preferred = ["Hu", "Agundez", "VULCAN", "velliet", "MOSESII", "stand", "agm2007"]

    def rank(n):
        return preferred.index(n) if n in preferred else 10_000 + discovered.index(n)

    net_order = sorted(discovered, key=rank)

    # map network name -> path
    net_path = {network_name_from_path(p): p for p in reprocess_rateconst_paths}

    # db:
    # key -> {
    #   "repr": reaction_text (first seen),
    #   "line": {net: ln or 0},
    #   "text": {net: reaction_text or ""},
    #   "tail": {net: tail_text or ""},  # from [50:]
    # }
    db = OrderedDict()

    RXN_W = 50

    for net in net_order:
        path = net_path[net]
        with open(path, "r") as f:
            for ln, raw in enumerate(f, 1):
                full = raw.rstrip("\n")

                # reaction is fixed-width in first 50 chars
                reaction_part = full[:RXN_W]
                reaction = reaction_part.strip()

                if not reaction:
                    continue
                if SEP_RXN not in reaction:
                    continue

                key = canonical_key(reaction)
                if key is None:
                    continue

                # tail is everything after col 50
                if "\t" in full:
                    tail = full.split("\t", 1)[1].rstrip()
                else:
                    tail = (full[RXN_W:] if len(full) > RXN_W else "").rstrip()
                if key not in db:
                    db[key] = {
                        "repr": reaction,
                        "line": {n: 0 for n in net_order},
                        "text": {n: "" for n in net_order},
                        "tail": {n: "" for n in net_order},
                    }

                if db[key]["line"][net] == 0:
                    db[key]["line"][net] = ln
                    db[key]["text"][net] = reaction
                    db[key]["tail"][net] = tail

    # -------------------------
    # write outputs
    # -------------------------
    out1 = f"uniq_reactions_rateconst_{tag}.dat"
    out2 = f"uniq_reactions_rateconst_metadata_{tag}.dat"
    out3 = f"uniq_reactions_rateconst_longtable_{tag}.dat"

    # file 1: unique reaction list
    with open(out1, "w") as f1:
        for entry in db.values():
            f1.write(entry["repr"] + "\n")

    # file 2: metadata + tail blocks
    COL_W = 10

    def fmt_rxn(s):
        s = (s or "").strip()
        return (s[:RXN_W]).ljust(RXN_W)

    def fmt_ln(n):
        try:
            n = int(n)
        except Exception:
            n = 0
        if n < 0:
            n = 0
        if n > 9999:
            n = 9999
        return f"{n:04d}".ljust(COL_W)

    header = [("Reaction".ljust(RXN_W))] + [net.ljust(COL_W) for net in net_order]

    tail_print_order = ["Hu", "Agundez", "VULCAN", "velliet", "MOSESII", "stand", "agm2007"]

    with open(out2, "w") as f2:
        f2.write("\t".join(header) + "\n")
        for entry in db.values():
            row = [fmt_rxn(entry["repr"])]
            row.extend(fmt_ln(entry["line"][net]) for net in net_order)
            f2.write("\t".join(row) + "\n")
            #print("entry line keys:", list(entry["line"].keys()))
            #print("tail_print_order:", tail_print_order)
            #break
            # Under each reaction, write tails in your requested order
            
            for net in tail_print_order:
                if net not in entry["line"]:
                    f2.write(f"{net}\t::\t!!something wrong with code!!\n")
                    continue

                if entry["line"][net] == 0:
                    f2.write(f"{net}\t::\t!!no reaction present!!\n")
                    continue

                tail = entry["tail"][net]
                if not tail.strip():
                    f2.write(f"{net}\t::\t!!no metadata present in org file!!\n")
                else:
                    f2.write(f"{net}\t::\t{tail}\n")

            f2.write("##########\n")  # spacer between reactions

    # file 3: optional long table (reaction text only, fixed 50 cols)
    if long_table == "Y":
        COL_W_LONG = 50

        def fmt_long(s):
            s = (s or "").strip()
            return s[:COL_W_LONG].ljust(COL_W_LONG)

        counts = {net: 0 for net in net_order}

        with open(out3, "w") as f3:
            f3.write("".join(fmt_long(net) for net in net_order) + "\n")

            for entry in db.values():
                row_parts = []
                for net in net_order:
                    txt = entry["text"][net]
                    if txt:
                        counts[net] += 1
                    row_parts.append(fmt_long(txt))
                f3.write("".join(row_parts) + "\n")

            f3.write("".join(fmt_long(str(counts[net])) for net in net_order) + "\n")

    print("done finding unique reactions + tails....")
    return {
        "tag": tag,
        "networks": net_order,
        "out_unique": out1,
        "out_metadata_with_tails": out2,
        "out_longtable": out3 if long_table == "Y" else None,
        "n_unique": len(db),
    }



#================================== UNIQUE SPECIES LIST & REACTANT NO SORTING ===============================

