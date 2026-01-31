import re
from pathlib import Path

#========================================= INPUTCHEM CURATION ACCORDING TO THE CHEMISTRY.DAT ======================
# INPUT :: agmfinal_netrxn.......
# OUTPUT :: agmfinal_input_curated....dat
# CHANGES :: just certain nomenclature changes to match chemistry.dat; + --> p, - --> a etc etc
# NUMBER OR REACTIONS SHOULD BE SAME IN BOTH
#==================================================================================================================

PLUS_SPLIT_EXACT = " + "  # IMPORTANT: separators are space+plus+space

# e^- bounded by spaces -> em
RE_EM = re.compile(r"(?<=\s)e\^-(?=\s)")

# HV bounded by spaces -> hnu
RE_HNU = re.compile(r"(?<=\s)HV(?=\s)")

# charge at END of token, like H^+, He^++, O2^- etc.
# left side: alphanumerics; right side: whitespace (token boundary in your format)
RE_CHARGE = re.compile(r"(?<=[A-Za-z0-9])\^(?P<sign>\+|-)(?P<rep>\+*|-*)(?=\s|\()")

def _out_name(input_path):
    name = Path(input_path).name
    tail = name[-20:]
    if tail.lower().endswith(".dat"):
        tail = tail[:-4]
    return f"agmfinal_input_curated_{tail}.dat"

def _replace_charges(s):
    def repl(m):
        sign = m.group("sign")
        rest = m.group("rep")  # extra + or - after the first one
        n = 1 + len(rest)
        if sign == "+":
            return "p" * n
        else:
            return "n" * n
    return RE_CHARGE.sub(repl, s)

def _compress_stoich(side_str):
    """
    side_str uses exact separator ' + '.
    Convert duplicates to 'n x TOKEN' (preserve first-seen order).
    """
    toks = [t.strip() for t in side_str.split(PLUS_SPLIT_EXACT) if t.strip()]

    counts = {}
    order = []
    for t in toks:
        if t not in counts:
            counts[t] = 0
            order.append(t)
        counts[t] += 1

    out = []
    for t in order:
        n = counts[t]
        if n > 1:
            out.append(f"{n} x {t}")
        else:
            out.append(t)

    return PLUS_SPLIT_EXACT.join(out)

def _curate_reaction_str(rxn):
    """
    rxn is ONLY the reaction text (the [:50] slice).
    Applies:
      - e^- -> em  (space-bounded)
      - HV  -> hnu (space-bounded)
      - ^+ ^++ -> p / pp, ^- ^-- -> a / aa  (end-of-token charge)
      - stoich compression on each side
    """
    s = rxn.rstrip("\n")

    # pad with spaces so the boundary rules work at ends too
    s = " " + s + " "

    s = RE_EM.sub("em", s)
    s = RE_HNU.sub("hnu", s)
    s = _replace_charges(s)

    # strip back padding
    s = s.strip()

    if "=" not in s:
        return s

    lhs, rhs = s.split("=", 1)
    lhs = lhs.strip()
    rhs = rhs.strip()

    # only split if it really uses " + " (your rule)
    lhs2 = _compress_stoich(lhs)
    rhs2 = _compress_stoich(rhs)

    return f"{lhs2} = {rhs2}"

def _process_line_style_file(input_path):
    """
    One reaction per line, reaction is line[:50]. Preserve the rest of the line unchanged.
    """
    input_path = Path(input_path)
    out_path = input_path.with_name(_out_name(input_path))

    with input_path.open("r", encoding="utf-8", errors="replace") as f, \
         out_path.open("w", encoding="utf-8") as w:
        for line in f:
            head = line[:50]
            tail = line[50:]  # includes tabs/columns/newline
            if "=" not in head:
                w.write(line)
                continue

            new_rxn = _curate_reaction_str(head)
            new_head = new_rxn[:50].ljust(50)
            w.write(new_head + tail)

    return out_path

def _process_kstack9_file(kpath):
    """
    One reaction per 9 lines. Reaction is in first line[:50].
    Keep all 9 lines; only modify the first line[:50].
    """
    kpath = Path(kpath)
    out_path = kpath.with_name(_out_name(kpath))

    with kpath.open("r", encoding="utf-8", errors="replace") as f, \
         out_path.open("w", encoding="utf-8") as w:

        block = []
        for line in f:
            block.append(line)
            if len(block) == 9:
                first = block[0]
                head = first[:50]
                tail = first[50:]

                if "=" in head:
                    new_rxn = _curate_reaction_str(head)
                    new_head = new_rxn[:50].ljust(50)
                    block[0] = new_head + tail

                for bline in block:
                    w.write(bline)
                block = []

        # if file ends with partial block, write it unchanged
        for bline in block:
            w.write(bline)

    return out_path

def agmfinal_master(input_path, kmode="N", kpath=None):
    """
    Master function.

    - Always processes input_path (one reaction per line, reaction in [:50]).
    - If kmode == 'Y', also processes kpath (one reaction per 9 lines).
    Output names:
      agmfinal_input_curated_{last14chars_of_input_filename}.dat
    Returns list of output paths.
    """
    outputs = []
    outputs.append(_process_line_style_file(input_path))

    if str(kmode).upper() == "Y":
        if kpath is None:
            raise ValueError("kmode='Y' but kpath was not provided")
        outputs.append(_process_kstack9_file(kpath))

    return outputs


file_path = "/Users/jb285991/Desktop/PHD/CODE/CHEMNETGEN/agmfinal_netrxn_addition_uniq_reactions_rateconst_H_HE_O.dat"
kfile_path = "/Users/jb285991/Desktop/PHD/CODE/CHEMNETGEN/agmfinal_netrxn_addition_uniq_reactions_rateconst_metadata_H_HE_O.dat"
# agmfinal_master("myfile.dat", kmode="N")
agmfinal_master(file_path, kmode="Y", kpath = kfile_path)
