import re
import os
from collections import Counter


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

def moses(path, allowed_elements, metadata="N"):
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
        species = re.sub(r'O\(\d+D\)', 'O', species)
        
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
            
            # First 45 characters as header
            if lineno == 1:
                kept.append("\n\n\n")
                kept.append(raw[:29])
                continue


            # Stop condition
            if stop_line_pattern in raw or lineno > 2115:
                break

            # Skip until "REACTIONS:" is found
            if not reactions_started:
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
    """
    Master function for the Agúndez chemical network.

    Parameters
    ----------
    path : str
        Path to Input_Agundez2025.dat
    allowed_elements : list
        Elements to keep, e.g. ['H','He','O']
    metadata : str
        'Y' to include metadata, 'N' otherwise
    """

    print("FILTERING AGUNDEZ CHEMICAL NETWORK....")

    allowed_upper = {el.upper() for el in allowed_elements}

    # -------------------------------------------------
    # Helpers
    # -------------------------------------------------

    def parse_formula(formula):
        """
        Parse a reduced chemical formula into elements.
        Example: C2H5O -> {'C':2,'H':5,'O':1}
        """
        tokens = re.findall(r'([A-Z][a-z]?)(\d*)', formula)
        comp = Counter()
        for el, n in tokens:
            comp[el.upper()] += int(n) if n else 1
        return comp

    def reduced_formula(species):
        """
        Reduce a species ONLY for elemental checking.
        Nothing here is written to output.
        """

        s = species.strip()

        # Always allowed helpers
        if s.upper() == "M":
            return None

        # Remove excited state after underscore (O_1D, CH2_1, etc)
        s = s.split("_")[0]

        # Remove any stray non-chemical symbols
        s = re.sub(r'[^A-Za-z0-9]', '', s)

        return s.strip()

    def species_allowed(species):
        reduced = reduced_formula(species)

        # M or helpers
        if reduced is None or not reduced:
            return True

        atoms = parse_formula(reduced)
        return set(atoms.keys()).issubset(allowed_upper)

    def extract_species(line):
        """
        Extract species from an Agúndez reaction line.
        Reactions use '=' as separator.
        """
        if "=" not in line:
            return []

        # Only reaction part, ignore kinetics/comments
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

    # -------------------------------------------------
    # Metadata control
    # -------------------------------------------------

    in_metadata = True
    discard_block = False

    # -------------------------------------------------
    # Read file
    # -------------------------------------------------

    with open(path, "r") as f:
        for line in f:
            raw = line.rstrip("\n")

            # Start of discarded reactions section (must be kept)
            if raw.startswith("! Reactions discarded because they involve species not included"):
                discard_block = True
                if metadata.upper() == "Y":
                    kept.append(raw)
                continue

            # References block starts → metadata resumes
            if raw.startswith("! # References"):
                in_metadata = True
                if metadata.upper() == "Y":
                    kept.append(raw)
                continue

            # Inside discarded reactions block → filter normally
            if discard_block and "=" in raw:
                species = extract_species(raw)
                if species and all(species_allowed(s) for s in species):
                    kept.append(raw)
                continue

            # Metadata lines
            if raw.startswith("!") and in_metadata:
                if metadata.upper() == "Y":
                    kept.append(raw)
                continue

            # Empty lines
            if not raw.strip():
                continue

            # Reaction filtering
            if "=" in raw:
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

    print(f"Filtered reactions written to {out_file}")
