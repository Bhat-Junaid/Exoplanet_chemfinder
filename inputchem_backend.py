from __future__ import annotations
from pathlib import Path
import re
from typing import Optional, Tuple

# Identifies the number of unique species in the final generated output file 
# filename: uniq_reactions_A_B_C.dat
def write_unique_species_list(path: str | Path, out_name: str = "unique_species_list.dat") -> tuple[int, Path]:
    """
    Extract unique species tokens from a reaction file and write them one per line.

    Output tokens are preserved (first-seen representative), e.g.:
      O2^+(X2Pig)
      HO4^+
      O2

    Uniqueness rules (canonical key):
    - Parentheses characters '(' and ')' are ignored (content kept):
        O(1D)  -> O1D
        O2^+(X2Pig) -> O2^+X2Pig
    - '^' is ignored for uniqueness (but '+' is kept if present):
        HO4^+ and HO4+ are the same key (if both appear)
    - Only '+' is considered part of the species if it is in the token text.
      (Separators are handled by splitting on spaced + and =)
    - Special case for tokens that are exactly two capital letters (optionally charged):
        AB == BA
        AB^+ == BA^+
        AB+ == BA+
    - Ignore token HV only when it is exactly 'HV'.

    Returns:
      (count_unique_species, output_path)
    """
    in_path = Path(path)
    if not in_path.exists():
        raise FileNotFoundError(in_path)

    out_path = in_path.with_name(out_name)
    
    # Split only on separators between species:
    # This assumes reactions are formatted with spaces around '+' and '=' like:
    #   A + B = C + D
    # It will NOT split the '+’ inside '^+' (no spaces there).
    split_re = re.compile(r"\s+(?:\+|=)\s+")

    # Clean tokens that might be quoted or have leading exclamation art
    lead_bang_re = re.compile(r"^!+\s*")
    trail_bang_re = re.compile(r"\s*!+$")

    def clean_token_for_output(tok: str) -> str:
        t = tok.strip()
        t = t.replace('"', '').strip()
        t = lead_bang_re.sub("", t).strip()
        t = trail_bang_re.sub("", t).strip()
        # collapse internal whitespace (shouldn’t exist inside species, but safe)
        t = re.sub(r"\s+", " ", t).strip()
        return t

    def canonical_key(token: str) -> str | None:
        """
        Build a canonical key used only for uniqueness.
        """
        t = token.strip()
        if not t:
            return None
        if t == "HV":
            return None

        # Ignore parentheses chars but keep content
        t = t.replace("(", "").replace(")", "")

        # For uniqueness: ignore '^' but keep '+'
        t = t.replace("^", "")

        # Special AB <-> BA equivalence ONLY for two-capital-letter tokens (optionally charged)
        # We interpret "two alphabets only" as exactly 2 uppercase letters in the base,
        # possibly followed by one or more '+'.
        m = re.fullmatch(r"([A-Z]{2})(\+*)", t)
        if m:
            base = m.group(1)
            charge = m.group(2)
            base = "".join(sorted(base))  # AB and BA map to same
            return base + charge

        return t

    # Map canonical_key -> representative original token (preserved)
    uniq: dict[str, str] = {}

    with in_path.open("r", encoding="utf-8", errors="replace") as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith("!"):
                continue

            # Split into raw tokens on spaced separators
            parts = split_re.split(line)

            for raw in parts:
                tok = clean_token_for_output(raw)
                if not tok:
                    continue
                if tok == "HV":
                    continue

                key = canonical_key(tok)
                if key is None:
                    continue

                # keep first seen representative (preserves e.g. O2^+(X2Pig))
                if key not in uniq:
                    uniq[key] = tok

    # Write in a stable order: sort by canonical key
    with out_path.open("w", encoding="utf-8") as f:
        for key in sorted(uniq.keys()):
            f.write(uniq[key] + "\n")

    return len(uniq), out_path


##########################
#########################
#########################


PLUS_SPLIT = re.compile(r"\s+\+\s+")

def sort_rxns_reactant(
    input_path: str | Path,
    kmode: str = "Y",
    kstack: int = 9,
    kpath: str | Path | None = None,
) -> Path:
    """
    Sort a .dat reaction file by number of reactants (ascending).

    Normal mode (kmode="N"): unchanged behavior
    ------------------------------------------
    Rules:
    - One reaction per line
    - Reactants are on the left of '=' or '=>'
    - Reactants are split by '+' only if surrounded by whitespace
    - 'HV' does NOT count as a reactant
    - Non-reaction lines are kept at the top, unchanged

    Output:
    - sorted_<inputfilename>.dat

    K-stack mode (kmode="Y"):
    -------------------------
    - Reads from kpath (NOT input_path)
    - Skips the FIRST line of kpath
    - Each reaction record is a block of kstack lines followed by a separator line '##########'
      (we detect blocks by the separator line; kstack is used as a sanity expectation)
    - Reaction used for sorting is block[0][:50]
    - After sorting blocks, writes the full block lines unchanged (including line1[50:] and the
      following lines), and also writes '##########' after each block.
    - Output filename is still: sorted_<kpath.name>.dat
    """

    input_path = Path(input_path)
    if input_path.suffix.lower() != ".dat":
        raise ValueError("Input file must be a .dat file")

    kmode = (kmode or "N").strip().upper()
    if kmode not in {"Y", "N"}:
        kmode = "N"

    def count_reactants(line: str) -> int | None:
        if "=>" in line:
            lhs = line.split("=>", 1)[0]
        elif "=" in line:
            lhs = line.split("=", 1)[0]
        else:
            return None

        parts = [p.strip() for p in PLUS_SPLIT.split(lhs) if p.strip()]
        return sum(1 for p in parts if p.upper() != "HV")

    # -------------------------
    # Normal mode: KEEP YOUR LOGIC
    # -------------------------
    if kmode == "N":
        output_path = input_path.with_name(f"sorted_{input_path.name}")

        headers: list[str] = []
        reactions: list[tuple[int, int, str]] = []

        with input_path.open("r", encoding="utf-8", errors="replace") as f:
            for i, raw in enumerate(f):
                line = raw.rstrip("\n")
                n = count_reactants(line)

                if n is None:
                    headers.append(line)
                else:
                    reactions.append((n, i, line))

        # stable sort: by reactant count, then original order
        reactions.sort(key=lambda x: (x[0], x[1]))

        with output_path.open("w", encoding="utf-8", newline="\n") as g:
            for line in headers:
                g.write(line + "\n")
            for _, _, line in reactions:
                g.write(line + "\n")

        return output_path

    # -------------------------
    # K-stack mode: sort blocks
    # -------------------------
    if kpath is None:
        raise ValueError("kmode='Y' requires kpath to be provided")

    kpath = Path(kpath)
    if kpath.suffix.lower() != ".dat":
        raise ValueError("kpath must be a .dat file")

    output_path = kpath.with_name(f"sorted_{kpath.name}")

    # Each block: (n_reactants, original_block_index, block_lines_list)
    blocks: list[tuple[int, int, list[str]]] = []
    header_lines: list[str] = []  # we skip first line entirely per your spec

    with kpath.open("r", encoding="utf-8", errors="replace") as f:
        # skip first line of the kpath file
        _ = next(f, None)

        current: list[str] = []
        block_idx = 0

        for raw in f:
            line = raw.rstrip("\n")

            if line.strip() == "##########":
                # end of one block
                if current:
                    # reaction is in first line, first 50 chars
                    rxn50 = current[0][:50].strip()
                    n = count_reactants(rxn50)

                    # If it doesn't look like a reaction, keep block as "header-ish" (top)
                    if n is None:
                        header_lines.extend(current)
                        header_lines.append("##########")
                    else:
                        blocks.append((n, block_idx, current + ["##########"]))
                        block_idx += 1

                current = []
                continue

            current.append(line)

        # handle file not ending with ##########
        if current:
            rxn50 = current[0][:50].strip()
            n = count_reactants(rxn50)
            if n is None:
                header_lines.extend(current)
            else:
                blocks.append((n, block_idx, current))
                block_idx += 1

    # stable sort blocks by (n_reactants, original_block_index)
    blocks.sort(key=lambda x: (x[0], x[1]))

    with output_path.open("w", encoding="utf-8", newline="\n") as g:
        # You asked to skip the first line; we do not re-write it.
        # If any non-reaction blocks exist (rare), they go first unchanged.
        for line in header_lines:
            g.write(line + "\n")

        for _, _, blines in blocks:
            for line in blines:
                g.write(line + "\n")

    return output_path

###########################
##########################
###########################


def rm_nonzero_agm(
    input_path: str | Path,
    clean_metadata: str = "Y",
    ksort: str = "Y",
    kpath: str | Path | None = None,
    kstack: int = 9,
) -> Tuple[Path, Optional[Path]]:
    """
    Two modes:

    A) ksort == 'Y' (default): stacked file mode (reaction every kstack lines)
       - reads from `kpath`
       - each block has `kstack` lines
       - the reaction line is the FIRST line of each block
       - keep the WHOLE block only if the LAST column of the reaction line == '0000'
       - output: netrxn_addition_<kpath.name>.dat
       - if no blocks pass, prints:
         "no 0000 entry in last column. Please check the files"
       - if clean_metadata == 'Y', ALSO write a "clean" file that contains ONLY the
         kept reaction lines (first line of each kept block), using the same meaning
         of clean_metadata as in mode B(just one change about keeping the whole line):
           * skip first line only if it starts with 'Reaction'
           * write exactly the whole line 
         Output name: netrxn_addition_<kpath.name with leading 'metadata_' removed>

    B) ksort == 'N': original one-line-reaction table mode
       1) Write a filtered file that removes rows where agm2007 != '0000'.
          Output: netrxn_addition_<inputfilename>

       2) If clean_metadata == 'Y' (default), also write a "clean" file built from
          the ORIGINAL INPUT, but ONLY for the rows that pass the agm2007 filter.
          - skips the first line ONLY if it starts with 'Reaction'
          - writes exactly 50 characters per kept line (truncate/pad)
          Output name: netrxn_addition_<inputfilename with leading 'metadata_' removed>
    """
    input_path = Path(input_path)

    yn = "Y" if str(clean_metadata).strip().upper().startswith("Y") else "N"
    kyn = "Y" if str(ksort).strip().upper().startswith("Y") else "N"

    def cleaned_basename(name: str) -> str:
        if name.startswith("metadata_"):
            return name[len("metadata_"):]
        return name.replace("metadata_", "", 1)

    # ----------------------------
    # Mode A: stacked blocks
    # ----------------------------
    if kyn == "Y":
        if kpath is None:
            raise ValueError("ksort='Y' requires kpath to be provided")

        kpath = Path(kpath)
        if kpath.suffix.lower() != ".dat":
            raise ValueError("kpath must be a .dat file")

        filtered_path = kpath.with_name(f"netrxn_addition_{kpath.name}")

        cleaned_path: Optional[Path] = None
        cleaned_fh = None

        if yn == "Y":
            cleaned_name = cleaned_basename(kpath.name)
            cleaned_path = kpath.with_name(f"netrxn_addition_{cleaned_name}")
            cleaned_fh = cleaned_path.open("w", encoding="utf-8", newline="\n")

        try:
            with kpath.open("r", encoding="utf-8", errors="replace") as f:
                all_lines = f.readlines()

            # header handling: skip first line only if it starts with "Reaction"
            has_header = bool(all_lines) and all_lines[0].lstrip().startswith("Reaction")
            if has_header:
                all_lines = all_lines[1:]

            kept_any = False

            with filtered_path.open("w", encoding="utf-8", newline="\n") as g:
                for i in range(0, len(all_lines), kstack):
                    block = all_lines[i : i + kstack]
                    if not block:
                        continue

                    # conservative: ignore incomplete final block
                    if len(block) < kstack:
                        continue

                    reaction_line = block[0]
                    stripped = reaction_line.strip()
                    if not stripped:
                        continue

                    cols = stripped.split()
                    if not cols:
                        continue

                    if cols[-1] == "0000":
                        kept_any = True
                        g.writelines(block)

                        # cleaned file in k-mode: write ONLY the reaction line
                        if cleaned_fh is not None:
                            s = reaction_line.rstrip("\n")
                            cleaned_fh.write(s + "\n")

            if not kept_any:
                print("no 0000 entry in last column. Please check the files")

        finally:
            if cleaned_fh is not None:
                cleaned_fh.close()

        return filtered_path, cleaned_path

    # ----------------------------
    # Mode B: original one-line table
    # ----------------------------
    filtered_path = input_path.with_name(f"netrxn_addition_{input_path.name}")

    cleaned_path: Optional[Path] = None
    cleaned_fh = None

    # detect header: skip first line only if it starts with 'Reaction'
    with input_path.open("r", encoding="utf-8", errors="replace") as f_probe:
        first_line = f_probe.readline()
    has_header = bool(first_line) and first_line.lstrip().startswith("Reaction")

    if yn == "Y":
        cleaned_name = cleaned_basename(input_path.name)
        cleaned_path = input_path.with_name(f"netrxn_addition_{cleaned_name}")
        cleaned_fh = cleaned_path.open("w", encoding="utf-8", newline="\n")

    try:
        with input_path.open("r", encoding="utf-8", errors="replace") as f, \
             filtered_path.open("w", encoding="utf-8", newline="\n") as g:

            for i, line in enumerate(f):
                stripped = line.strip()

                # always keep empty lines (in filtered)
                if not stripped:
                    g.write(line)
                    continue

                cols = stripped.split()

                # header or malformed line → keep in filtered
                if cols[-1].isalpha() or len(cols[-1]) != 4:
                    g.write(line)
                    continue

                # keep only agm2007 == 0000
                if cols[-1] == "0000":
                    g.write(line)

                    # cleaned file: skip first line only if it is a header line
                    if cleaned_fh is not None and not (has_header and i == 0):
                        s = line.rstrip("\n")
                        cleaned_fh.write(s[:50].ljust(50) + "\n")

    finally:
        if cleaned_fh is not None:
            cleaned_fh.close()

    return filtered_path, cleaned_path


""""
def rm_nonzero_agm1(input_path: str | Path, clean_metadata: str = "Y") -> Tuple[Path, Optional[Path]]:

    1) Write a filtered file that removes rows where agm2007 != '0000'.
       Output: netrxn_addition_<inputfilename>

    2) If clean_metadata == 'Y' (default), also write a "clean" file built from
       the ORIGINAL INPUT, but ONLY for the rows that pass the agm2007 filter.
       - skips the first line of the input
       - writes exactly 50 characters per kept line (truncate/pad)
       Output name: netrxn_addition_<inputfilename with leading 'metadata_' removed>
                   e.g. metadata_H_HE_O.dat -> netrxn_addition_H_HE_O.dat

    input_path = Path(input_path)
    yn = "Y" if str(clean_metadata).strip().upper().startswith("Y") else "N"

    filtered_path = input_path.with_name(f"netrxn_addition_{input_path.name}")

    cleaned_path: Optional[Path] = None
    cleaned_fh = None

    if yn == "Y":
        cleaned_name = input_path.name
        if cleaned_name.startswith("metadata_"):
            cleaned_name = cleaned_name[len("metadata_"):]
        else:
            cleaned_name = cleaned_name.replace("metadata_", "", 1)

        cleaned_path = input_path.with_name(f"netrxn_addition_{cleaned_name}")
        cleaned_fh = cleaned_path.open("w", encoding="utf-8", newline="\n")

    try:
        with input_path.open("r", encoding="utf-8", errors="replace") as f, \
             filtered_path.open("w", encoding="utf-8", newline="\n") as g:

            for i, line in enumerate(f):
                stripped = line.strip()

                # always keep empty lines (in filtered); for cleaned, write 50 cols if kept
                if not stripped:
                    g.write(line)
                    continue

                cols = stripped.split()

                # header or malformed line → keep in filtered, but cleaned skips first line anyway
                if cols[-1].isalpha() or len(cols[-1]) != 4:
                    g.write(line)
                    continue

                # keep only agm2007 == 0000
                if cols[-1] == "0000":
                    g.write(line)

                    # cleaned file: from ORIGINAL INPUT, skip first line of the file
                    if cleaned_fh is not None and i != 0:
                        s = line.rstrip("\n")
                        cleaned_fh.write(s[:49].ljust(49) + "\n")

    finally:
        if cleaned_fh is not None:
            cleaned_fh.close()

    return filtered_path, cleaned_path


"""