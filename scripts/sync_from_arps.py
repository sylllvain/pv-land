#!/usr/bin/env python3
"""Synchronise the PV-LAND offline model with the ARPS PV-LAND sources.

ARPS is the reference:
  - core/ and arps_lib/ are exact copies of ARPS files (or of complete
    subroutines extracted from ARPS files);
  - offline/pvland_sfcphy_offline.f90 is generated from the PV-LAND block of
    ARPS sfcphy3d.f90 (sfcphy == 5) inserted into the template
    scripts/pvland_sfcphy_offline.tmpl.f90, with a few lines tagged [OFFLINE].
These files must never be edited by hand here: edit ARPS (or the template),
then re-run this script.

Usage (from the PV-LAND directory):
    python3 scripts/sync_from_arps.py            # copy/generate from ARPS
    python3 scripts/sync_from_arps.py --check    # only compare, exit 1 if different

The ARPS directory defaults to ../arps5.3.4_PVLAND (override with --arps DIR).
"""
import argparse
import difflib
import os
import re
import sys

# Whole files copied verbatim from <ARPS>/src/arps
CORE_FILES = [
    "module_pvland.f90",
    "pvland.f90",
    "pvland_profil.f90",
    "pvland_view_factors.f90",
    "pvland_solar_radiation.f90",
    "pvland_radiosity_SW.f90",
    "pvland_radiosity_LW.f90",
    "pvland_pv_thermal.f90",
]
ARPS_LIB_FILES = [
    "soilebm3d.f90",
    "thermolib3d.f90",
]
# Complete subroutines extracted from larger ARPS files:
#   output file in arps_lib/  ->  (ARPS source file, [subroutines, in order])
ARPS_LIB_EXTRACTS = {
    "sfcphy3d_extract.f90": ("sfcphy3d.f90", [
        "sfcflx", "dragcoef", "cuc", "cptc",
        "cucwtr", "cptcwtr", "cuneuwtr", "cptneuwtr"]),
    "tmix3d_extract.f90": ("tmix3d.f90", ["tridiag2"]),
}
# Include files copied verbatim from <ARPS>/include
INCLUDE_FILES = [
    "globcst.inc",
    "phycst.inc",
    "sfcphycst.inc",
    "soilcst.inc",
    "grid.inc",
]

HERE = os.path.dirname(os.path.abspath(__file__))
STEP_TEMPLATE = os.path.join(HERE, "pvland_sfcphy_offline.tmpl.f90")
STEP_MARKER = "!@@ARPS_SFCPHY5_BLOCK@@\n"


def read(path):
    return open(path, encoding="latin-1").read()


def extract_subroutines(path, names):
    """Return the text of the given subroutines (comment banner included)."""
    lines = read(path).splitlines(keepends=True)
    starts = {}
    for i, line in enumerate(lines):
        m = re.match(r"\s*SUBROUTINE\s+(\w+)", line, re.I)
        if m:
            starts.setdefault(m.group(1).lower(), i)
    out = []
    for name in names:
        if name not in starts:
            sys.exit(f"ERROR: SUBROUTINE {name} not found in {path}")
        i0 = starts[name]
        # include the comment banner just above the SUBROUTINE statement
        while i0 > 0 and (lines[i0 - 1].lstrip().startswith("!")
                          or not lines[i0 - 1].strip()):
            i0 -= 1
        i1 = starts[name]
        while not re.match(r"\s*END\s+SUBROUTINE", lines[i1], re.I):
            i1 += 1
        out.append(f"!==== Extracted from ARPS {os.path.basename(path)}"
                   f" : SUBROUTINE {name} (do not edit) ====\n")
        out.extend(lines[i0:i1 + 1])
        out.append("\n")
    return "".join(out)


def comment_out(text):
    """Comment out every non-empty line (keeps the original ARPS code visible)."""
    return "".join(("!" + line[1:] if line[:1] == " " else "!" + line) if line.strip() else line
                   for line in text.splitlines(keepends=True))


def make_offline_step(sfcphy_path):
    """Offline driver = ARPS sfcphy==5 block + lines tagged [OFFLINE]."""
    lines = read(sfcphy_path).splitlines(keepends=True)
    try:
        i0 = next(i for i, l in enumerate(lines)
                  if re.match(r"\s*ELSE IF \( sfcphy == 5 \) THEN", l)) + 1
        i1 = next(i for i, l in enumerate(lines)
                  if i > i0 and "DEALLOCATE(lambda_deg,phi_deg)" in l)
    except StopIteration:
        sys.exit(f"ERROR: PV-LAND block (sfcphy == 5) not found in {sfcphy_path}")
    body = "".join(lines[i0:i1 + 1])

    def substitute(pattern, new_lines, keep_old=True):
        nonlocal body
        found = re.findall(pattern, body, re.S)
        if len(found) != 1:
            sys.exit(f"ERROR: ARPS PV-LAND block changed, cannot find exactly once:\n  {pattern}\n"
                     f"  -> update make_offline_step() in {__file__}")
        old = found[0]
        body = body.replace(old, (comment_out(old) if keep_old else "") + new_lines)

    # 1. direct/diffuse split of SW -> measured diffuse SW
    substitute(r"( *! split solar radiation between direct and diffuse\n"
               r" *CALL split_sw\(.*?1367\.\)\n)",
               "      S_diffuse(i,j) = S_diffuse_in(i,j)                 ! [OFFLINE] measured diffuse SW\n"
               "      S_direct(i,j)  = radsw(i,j) - S_diffuse(i,j)       ! [OFFLINE]\n")
    # 2. solar angles in the park frame -> measured lambda, phi (rad)
    substitute(r"( *call compute_lambda_phi\(szenith_deg[^\n]*\n)",
               "      lambda_deg(i,j) = lambda_in(i,j) * 180. / pi       ! [OFFLINE] measured, park frame\n"
               "      phi_deg(i,j)    = phi_in(i,j)    * 180. / pi       ! [OFFLINE] measured, park frame\n")
    # 3. Obukhov length from the surface fluxes -> measured Obukhov length
    substitute(r"(\s*if \(ptsflx\(i,j\)\.ne\.0\.\) then\n.*?endif\n)",
               "      L_obk(i,j) = L_obk_in(i,j)                         ! [OFFLINE] measured Obukhov length\n")
    # 4. save local diagnostics before they are deallocated
    substitute(r"(  DEALLOCATE\(cdm_u,cdm_v,z0_PV\)\n)",
               "  uref_diag    = uref        ! [OFFLINE DIAG]\n"
               "  cdh_can_diag = cdh_can     ! [OFFLINE DIAG]\n"
               "  cdh_PVa_diag = cdh_PVa     ! [OFFLINE DIAG]\n\n"
               "  DEALLOCATE(cdm_u,cdm_v,z0_PV)\n", keep_old=False)

    template = read(STEP_TEMPLATE)
    if template.count(STEP_MARKER) != 1:
        sys.exit(f"ERROR: marker {STEP_MARKER.strip()} not found in {STEP_TEMPLATE}")
    return template.replace(STEP_MARKER, body)


def targets(arps):
    """Yield (destination, expected text) for every synchronised file."""
    src = os.path.join(arps, "src", "arps")
    inc = os.path.join(arps, "include")
    for f in CORE_FILES:
        yield os.path.join("core", f), read(os.path.join(src, f))
    for f in ARPS_LIB_FILES:
        yield os.path.join("arps_lib", f), read(os.path.join(src, f))
    for out, (f, names) in ARPS_LIB_EXTRACTS.items():
        yield os.path.join("arps_lib", out), extract_subroutines(os.path.join(src, f), names)
    for f in INCLUDE_FILES:
        yield os.path.join("arps_lib", "include", f), read(os.path.join(inc, f))
    yield (os.path.join("offline", "pvland_sfcphy_offline.f90"),
           make_offline_step(os.path.join(src, "sfcphy3d.f90")))


def main():
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--arps", default=os.path.join("..", "arps5.3.4_PVLAND"),
                   help="ARPS PV-LAND directory (default: ../arps5.3.4_PVLAND)")
    p.add_argument("--check", action="store_true",
                   help="compare only, do not write; exit 1 if files differ")
    args = p.parse_args()
    if not os.path.isdir(os.path.join(args.arps, "src", "arps")):
        sys.exit(f"ERROR: ARPS directory not found: {args.arps}")

    ndiff = 0
    for dest, text in targets(args.arps):
        current = read(dest) if os.path.exists(dest) else None
        if current == text:
            continue
        ndiff += 1
        if args.check:
            print(f"DIFFERENT: {dest}")
            diff = difflib.unified_diff(
                (current or "").splitlines(), text.splitlines(),
                f"{dest} (offline)", f"{dest} (from ARPS)", lineterm="", n=1)
            for k, line in enumerate(diff):
                if k >= 40:
                    print("   ...")
                    break
                print("   " + line)
        else:
            os.makedirs(os.path.dirname(dest), exist_ok=True)
            with open(dest, "w", encoding="latin-1") as fh:
                fh.write(text)
            print(f"updated: {dest}")

    if args.check:
        if ndiff:
            print(f"{ndiff} file(s) differ from ARPS ({args.arps}).")
            sys.exit(1)
        print(f"OK: core/, arps_lib/ and offline/pvland_sfcphy_offline.f90 match ARPS ({args.arps}).")
    elif ndiff == 0:
        print("Nothing to update: already identical to ARPS.")


if __name__ == "__main__":
    main()
