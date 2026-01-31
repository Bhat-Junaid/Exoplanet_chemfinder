from pathlib import Path
import re

PLUS_SPLIT = re.compile(r"\s*\+\s*")

def strip_parens(tok: str) -> str:
    return tok.replace("(", "").replace(")", "")

def normalize_token_for_matching(tok: str) -> str:
    """
    Matching-only normalization:
    - remove '(' and ')'
    - if the resulting token is ONLY uppercase letters (A–Z), sort letters (OH == HO)
    - otherwise leave unchanged (digits, charges, excited-state text, etc.)
    """
    t = strip_parens(tok.strip())
    if re.fullmatch(r"[A-Z]+", t):
        t = "".join(sorted(t))
    return t


def parse_reaction_line(line):
    raw = line.strip()
    if not raw or raw.startswith("#"):
        return None
    if "=" not in raw:
        return None

    lhs, rhs = raw.split("=", 1)

    lhs_tokens = [normalize_token_for_matching(t) for t in PLUS_SPLIT.split(lhs) if t.strip()]
    rhs_tokens = [normalize_token_for_matching(t) for t in PLUS_SPLIT.split(rhs) if t.strip()]

    lhs_key = tuple(sorted(lhs_tokens))
    rhs_key = tuple(sorted(rhs_tokens))

    # Direction-sensitive: A=B is NOT the same as B=A
    return (lhs_key, rhs_key)


def load_reactions(path: str):
    """
    Returns:
    - keys: set of canonical keys
    - key_to_line: representative original line for each key (first seen)
    """
    keys = set()
    key_to_line = {}

    with open(path, "r", encoding="utf-8", errors="replace") as f:
        for line in f:
            k = parse_reaction_line(line)
            if k is None:
                continue
            keys.add(k)
            if k not in key_to_line:
                key_to_line[k] = line.rstrip("\n")

    return keys, key_to_line



def diff_reactions(file1, file2, out_path=None):
    """
    Master function.

    Inputs
    - file1, file2: paths (str or Path)
    - out_path: optional output path (str or Path). If None -> diff_<stem1>__<stem2>.dat

    Output
    - Writes the symmetric-difference reactions to out_path
    - Returns a dict with counts + output path
    """
    file1 = Path(file1)
    file2 = Path(file2)

    keys1, map1 = load_reactions(str(file1))
    keys2, map2 = load_reactions(str(file2))

    only_in_1 = keys1 - keys2
    only_in_2 = keys2 - keys1
    diff_keys = sorted(only_in_1 | only_in_2)

    out_lines = []
    for k in diff_keys:
        if k in map1:
            out_lines.append(map1[k])
        else:
            out_lines.append(map2[k])

    out_lines.sort()

    if out_path is None:
        out_path = file1.with_name(f"diff_{file1.stem}__{file2.stem}.dat")
    else:
        out_path = Path(out_path)

    with out_path.open("w", encoding="utf-8") as w:
        for line in out_lines:
            w.write(line + "\n")

    summary = {
        "file1": str(file1),
        "file2": str(file2),
        "file1_unique": len(keys1),
        "file2_unique": len(keys2),
        "only_in_file1": len(only_in_1),
        "only_in_file2": len(only_in_2),
        "symmetric_diff": len(diff_keys),
        "output": str(out_path),
    }

    print(
        f"file1 unique: {summary['file1_unique']}\n"
        f"file2 unique: {summary['file2_unique']}\n"
        f"only in file1: {summary['only_in_file1']}\n"
        f"only in file2: {summary['only_in_file2']}\n"
        f"symmetric diff (not in both): {summary['symmetric_diff']}\n"
        f"wrote: {summary['output']}"
    )

    return summary

diff_reactions("/Users/jb285991/Desktop/PHD/CODE/CHEMNETGEN/uniq_reactions_rateconst_H_HE_O.dat",
               "/Users/jb285991/Desktop/PHD/CODE/CHEMNETGEN/uniq_reactions_H_HE_O.dat",
               "/Users/jb285991/Desktop/PHD/CODE/CHEMNETGEN/reaction_differences.dat"
               )