import re
from pathlib import Path

# numeric metadata columns look like 0000, 0184, ...
NUM4 = re.compile(r"^\d{4}$")

# split reaction around + and arrow
PLUS_SPLIT = re.compile(r"\s*\+\s*")

def find_numeric_start(parts):
    for i, tok in enumerate(parts):
        if NUM4.match(tok):
            return i
    return None

def normalize_arrow(rxn: str) -> str:
    # support '=' and '->' and '=>'
    rxn = rxn.replace("=>", "->").replace("=", "->")
    return " ".join(rxn.split())

def species_to_latex(sp: str) -> str:
    """
    Convert a species token like:
      H3^+   -> H_3^{+}
      He^++  -> He^{2+}
      e^-    -> e^{-}
      O1D    -> O(^1D)  (common journal style)
      HV     -> h\\nu
      M      -> M
    Keeps things conservative (doesn't guess weird formats).
    """
    sp = sp.strip()

    # special common cases
    if sp.upper() == "HV":
        return r"h\nu"
    if sp == "M":
        return "M"

    # electron
    if sp.startswith("e^"):
        ch = sp.split("^", 1)[1]
        return rf"e^{{{ch}}}"

    # handle charge written as ^+, ^-, ^++, ^--, ^+2, ^-2, etc.
    charge = None
    if "^" in sp:
        base, ch = sp.split("^", 1)
        sp = base
        charge = ch.strip()

    # excited oxygen shorthand like O1D, O1S (common in atmospheric networks)
    # convert to O(^1D) / O(^1S) if it matches that pattern
    m_exc = re.fullmatch(r"([A-Za-z]+)(\d)([A-Za-z])", sp)
    if m_exc:
        el, mult, state = m_exc.groups()
        base_tex = rf"{el}(^{{{mult}{state}}})"
    else:
        # standard chemical formula: element symbols + optional digits
        # e.g. H3O -> H_3O, H2O2 -> H_2O_2, He -> He
        base_tex = re.sub(r"(\d+)", r"_\1", sp)

    # now apply charge
    if charge is None:
        return base_tex

    # normalize repeated signs: '++' -> '2+'
    if set(charge) <= {"+", "-"} and len(charge) > 1:
        charge = f"{len(charge)}{charge[0]}"

    return rf"{base_tex}^{{{charge}}}"

def reaction_to_math_latex(rxn: str) -> str:
    rxn = normalize_arrow(rxn)

    if "->" not in rxn:
        # fallback: just render whole thing as \mathrm{...}
        return rf"$\mathrm{{{rxn}}}$"

    left, right = [s.strip() for s in rxn.split("->", 1)]

    def side_to_tex(side: str) -> str:
        terms = PLUS_SPLIT.split(side) if side else []
        terms_tex = []
        for t in terms:
            t = t.strip()
            if not t:
                continue
            terms_tex.append(species_to_latex(t))
        return " + ".join(terms_tex)

    L = side_to_tex(left)
    R = side_to_tex(right)

    return rf"$\mathrm{{{L} \rightarrow {R}}}$"

def parse_line(line: str):
    line = line.rstrip("\n")
    if not line.strip():
        return None
    parts = line.split()  # tabs/spaces OK
    idx = find_numeric_start(parts)
    if idx is None:
        return None
    rxn = " ".join(parts[:idx]).strip()
    cols = parts[idx:]
    return rxn, cols

# --------- paths ----------
infile = Path("/Users/jb285991/Desktop/PHD/CODE/CHEMNETGEN/Uniq_reactions_metadata_H_HE_O.dat")  # <-- change to your file path/name
outfile = Path("reactions_longtable.tex")

rows = []
max_cols = 0

with infile.open("r", encoding="utf-8") as f:
    for line in f:
        parsed = parse_line(line)
        if not parsed:
            continue
        rxn, cols = parsed
        rows.append((rxn, cols))
        max_cols = max(max_cols, len(cols))

# pad numeric columns so LaTeX table is rectangular
rows = [(rxn, cols + [""] * (max_cols - len(cols))) for rxn, cols in rows]

# build LaTeX longtable
col_spec = "l" + "r" * max_cols

tex = []
tex.append(r"\begin{longtable}{" + col_spec + r"}")
tex.append(r"\hline")
tex.append(
    r"\textbf{Reaction} & " +
    " & ".join([rf"\textbf{{C{i+1}}}" for i in range(max_cols)]) +
    r" \\"
)
tex.append(r"\hline")
tex.append(r"\endfirsthead")

tex.append(r"\hline")
tex.append(
    r"\textbf{Reaction (cont.)} & " +
    " & ".join([rf"\textbf{{C{i+1}}}" for i in range(max_cols)]) +
    r" \\"
)
tex.append(r"\hline")
tex.append(r"\endhead")

tex.append(r"\hline")
tex.append(r"\endfoot")

for rxn, cols in rows:
    rxn_tex = reaction_to_math_latex(rxn)
    tex.append(rxn_tex + " & " + " & ".join(cols) + r" \\")

tex.append(r"\hline")
tex.append(r"\end{longtable}")

outfile.write_text("\n".join(tex) + "\n", encoding="utf-8")
print(f"Wrote: {outfile}")
