from __future__ import annotations
from pathlib import Path
import re


def agmchemnet(path: str | Path, skip_lines: int = 40, rate_const: str = "Y") -> Path:
    """
    Master function for the Antonio/AGM chemistry file reprocessing.

    Pipeline
    1) Skip first `skip_lines` lines.
    2) Detect reaction blocks via lines containing 'NEW REACTION' (case-insensitive).
    3) For each marker:
       - extract reaction line exactly 3 lines below it
       - if rate_const == 'Y', also extract rate-constant line exactly 6 lines below it
    4) Clean and normalize reaction line:
       - remove quotes and exclamation art
       - '->' -> '='
       - whole-word 'hnu' (any case) -> 'HV'
       - 'O(3p)' -> 'O(3P)'
       - then 'O(3P)' -> 'O'
       - 'O2(ground)' -> 'O2'
       - 'em' -> 'e^-'
       - every literal 'p' (case-sensitive) -> '^+'
       - expand patterns like 'N x TOKEN' into 'TOKEN + TOKEN + ...' (N times),
         where ' x ' is exactly space-x-space.

    Output
    - If rate_const == 'N':
        writes one reaction per line to: reprocess_agm2007_output_CO2.dat
    - If rate_const == 'Y':
        writes: reaction (50 chars, left-justified) + rate-constant line (as-is, stripped)
        to: reprocess_rate_const_agm2007_output_CO2.dat
    """
    print("REPROCESSING ANTONIO/AGM CHEMISTRY FILE...")

    in_path = Path(path)
    if in_path.suffix.lower() != ".dat":
        raise ValueError("Expected a .dat input file")

    rc = rate_const.strip().upper()
    if rc not in {"Y", "N"}:
        raise ValueError("rate_const must be 'Y' or 'N'")

    if rc == "Y":
        out_path = in_path.with_name("reprocess_rate_const_agm2007_output_CO2.dat")
    else:
        out_path = in_path.with_name("reprocess_agm2007_output_CO2.dat")

    # -------------------------------------------------
    # Patterns / helpers
    # -------------------------------------------------
    new_reaction_re = re.compile(r"\bNEW\s+REACTION\b", re.IGNORECASE)
    hnu_re = re.compile(r"\bhnu\b", re.IGNORECASE)

    # N x TOKEN where x has spaces around it; TOKEN = non-space run
    mult_re = re.compile(r"\b(\d+)\s+x\s+(\S+)")

    def expand_multipliers(expr: str) -> str:
        s = expr
        while True:
            m = mult_re.search(s)
            if not m:
                break
            n = int(m.group(1))
            tok = m.group(2)
            repl = " + ".join([tok] * n) if n > 0 else ""
            s = s[:m.start()] + repl + s[m.end():]
        return s

    def clean_reaction_line(line: str) -> str:
        s = line.replace('"', '').strip()

        # remove leading runs of exclamation marks and whitespace
        s = re.sub(r"^!+\s*", "", s).strip()

        # remove trailing exclamation art
        s = s.strip("!").strip()

        # base conversions
        s = s.replace("->", "=")
        s = hnu_re.sub("HV", s)

        # requested string normalizations (order matters)
        s = s.replace("O(3p)", "O(3P)")
        s = s.replace("O(3P)", "O")
        s = s.replace("O2(ground)", "O2")
        s = s.replace("em", "e^-")

        # convert every literal 'p' (case-sensitive) to '^+'
        s = s.replace("p", "^+")

        # normalize whitespace before expanding multipliers
        s = re.sub(r"\s+", " ", s).strip()

        # expand 'N x TOKEN'
        s = expand_multipliers(s)

        # final whitespace cleanup
        s = re.sub(r"\s+\+\s+", " + ", s)
        s = re.sub(r"\s*=\s*", " = ", s).strip()

        return s

    def clean_rate_line(line: str) -> str:
        # keep it basically as-is, just strip newline and outer spaces
        return line.rstrip("\n").strip()

    # -------------------------------------------------
    # Read + process
    # -------------------------------------------------
    with in_path.open("r", encoding="utf-8", errors="replace") as f:
        lines = f.readlines()

    lines = lines[min(skip_lines, len(lines)) :]

    out_lines: list[str] = []
    i = 0
    nlines = len(lines)

    while i < nlines:
        if new_reaction_re.search(lines[i]):
            rxn_idx = i + 3
            k_idx = i + 6  # rate constant line

            if rxn_idx < nlines:
                raw_rxn = lines[rxn_idx]
                rxn = clean_reaction_line(raw_rxn)

                if rxn and "=" in rxn:
                    if rc == "Y":
                        kline = clean_rate_line(lines[k_idx]) if k_idx < nlines else ""
                        out_lines.append(rxn.ljust(50) +"\t" + kline + "\n")
                    else:
                        out_lines.append(rxn + "\n")

            i = rxn_idx + 1
            continue

        i += 1

    # -------------------------------------------------
    # Write output
    # -------------------------------------------------
    with out_path.open("w", encoding="utf-8") as f:
        f.writelines(out_lines)

    print(f"DONE. Wrote {len(out_lines)} reactions to: {out_path}")
    return out_path

agmchemnet("/Users/jb285991/Desktop/PHD/CODE/CHEMNETGEN/input_agm2007.dat")
