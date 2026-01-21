
from pathlib import Path
import re

from collections import Counter

PLUS_SPLIT = re.compile(r" \+ ")

"""
For the verification of number of species;
verification of unique_species_list.dat
It takes the uniq_reactions_agm.dat and splits each 
species of the equation into a single line and then 
finds the repeated lines or unique lines. Thus giving 
unique species.

*********
CAUTION:
********
This is a preliminary code. It includes HV as species 
which is not the case as it is just a proxy for light
This code should only be used as a sanity check to 
unique_species_list.dat not the primary source to generate 
such kind of species list.

"""


# Two capital letters, optional charge with optional caret: AB, BA, AB+, AB^+, AB-, AB^-
_TWO_LETTER_CHARGED_RE = re.compile(r"^([A-Z])([A-Z])(?:\^?([+-]))?$")

def reactions_to_tokens(
    input_path: str | Path,
    output_path: str | Path,
):
    """
    Read a .dat reaction file and write one species token per line.

    Rules:
    - One reaction per line
    - Reactants + products are split by '=' or '=>'
    - Species are split by ' + ' (exactly one space, plus, one space)
    - Only '(' and ')' characters are removed (O(1D) -> O1D)
    - FINAL step: if token is exactly two capital letters (optionally with +/- charge),
      reorder letters alphabetically so AB == BA, AB^+ == BA^+, AB- == BA-, etc.
    """

    input_path = Path(input_path)
    output_path = Path(output_path)

    if input_path.suffix.lower() != ".dat":
        raise ValueError("Input file must be a .dat file")

    def strip_paren_chars(token: str) -> str:
        return token.replace("(", "").replace(")", "")

    def normalize_two_letter_order(token: str) -> str:
        """
        Apply ONLY to tokens like AB, BA, AB^+, BA+, AB^-, etc.
        Keep the charge sign as-is, keep whether caret exists as-is.
        """
        m = _TWO_LETTER_CHARGED_RE.fullmatch(token)
        if not m:
            return token

        a, b, charge = m.group(1), m.group(2), m.group(3)

        # canonical order: alphabetical
        base = "".join(sorted([a, b]))

        if charge is None:
            return base

        # preserve caret presence if it was written with caret in the original token
        has_caret = "^" in token
        return f"{base}{'^' if has_caret else ''}{charge}"

    def tokenize_side(side: str) -> list[str]:
        side = side.strip()
        if not side:
            return []

        parts = [p.strip() for p in side.split(" + ") if p.strip()]
        out = []
        for p in parts:
            p = strip_paren_chars(p)
            out.append(p)
        return out

    tokens: list[str] = []

    with input_path.open("r", encoding="utf-8", errors="replace") as f:
        for line in f:
            line = line.strip()
            if not line:
                continue

            if "=>" in line:
                lhs, rhs = line.split("=>", 1)
            elif "=" in line:
                lhs, rhs = line.split("=", 1)
            else:
                continue

            tokens.extend(tokenize_side(lhs))
            tokens.extend(tokenize_side(rhs))

    # FINAL pass: apply AB==BA rule after tokens exist as lines
    tokens = [normalize_two_letter_order(tok) for tok in tokens]

    output_path.parent.mkdir(parents=True, exist_ok=True)
    with output_path.open("w", encoding="utf-8") as f:
        for tok in tokens:
            f.write(tok + "\n")

reactions_to_tokens(
    "/Users/jb285991/Desktop/PHD/CODE/CHEMNETGEN/Uniq_reactions_agm_H_HE_O.dat",
    "tokens_out.dat"
)



filename = "/Users/jb285991/Desktop/PHD/CODE/CHEMNETGEN/tokens_out.dat"

with open(filename, "r") as f:
    lines = [line.rstrip("\n") for line in f]

# count occurrences
counts = Counter(lines)

# total stats
total_lines = len(lines)
unique_lines = len(counts)

# filter only repeated lines
repeated = {line: count for line, count in counts.items() if count >= 1}

print(f"Total lines        : {total_lines}")
print(f"Unique lines       : {unique_lines}")
print(f"Repeated line types: {len(repeated)}\n")

if not repeated:
    print("No repeated lines found.")
else:
    print("Repeated lines:\n")
    for line, count in sorted(repeated.items()):
        print(f"{count} times  |  {line}")


