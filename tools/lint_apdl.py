"""
Static checks for the APDL macros (no ANSYS needed).

  * *DO/*ENDDO and *IF..THEN/*ENDIF balance (with line numbers)
  * command line longer than 640 characters (APDL limit)
  * parameter names: <= 32 characters, start with a letter
  * parameter names that clash with APDL get-functions or DOF labels
    (NX, NY, UY, TEMP, NODE, NSEL ...) or with a component (CM) name
  * across macros: a component name that is a parameter in another macro
    (RESUME restores parameters)
  * MPDATA / TBDATA with more than 6 values, array assignment with more
    than 10 values (per-command limits)
  * *VWRITE format lines longer than 80 characters

Usage:  python3 tools/lint_apdl.py ansys/*.mac
"""
import re
import sys

RESERVED = {"NX", "NY", "NZ", "UX", "UY", "UZ", "TEMP", "NODE", "NSEL", "ESEL", "KSEL",
            "LSEL", "ASEL", "VSEL", "CENTRX", "CENTRY", "ROTX", "ROTY", "ROTZ", "PRES",
            "VOLT", "ALL", "STAT", "DEFA", "LOC", "ANG", "KP", "LINE", "AREA", "VOLU",
            "ELEM", "SIN", "COS", "TAN", "EXP", "LOG", "ABS", "SQRT", "NINT", "MOD", "MAX", "MIN"}


ALL_PARAMS, ALL_COMPS = {}, {}


def lint(path):
    issues, stack = [], []
    params, comps = {}, {}
    lines = open(path).read().splitlines()
    for no, raw in enumerate(lines, 1):
        code = raw.split("!")[0].rstrip()
        s = code.strip()
        if not s:
            continue
        if len(raw) > 640:
            issues.append((no, "line longer than 640 characters"))
        up = s.upper()
        f = [x.strip() for x in s.split(",")]
        cmd = f[0].upper()
        if cmd == "*DO":
            stack.append(("*DO", no))
        elif cmd == "*ENDDO":
            if not stack or stack[-1][0] != "*DO":
                issues.append((no, "*ENDDO without *DO"))
            else:
                stack.pop()
        elif cmd == "*IF" and f[-1].upper() == "THEN":
            stack.append(("*IF", no))
        elif cmd == "*ENDIF":
            if not stack or stack[-1][0] != "*IF":
                issues.append((no, "*ENDIF without *IF"))
            else:
                stack.pop()
        m = re.match(r"^([A-Za-z_][A-Za-z0-9_]*)(\([^)]*\))?\s*=\s*(.+)$", s)
        if m and not up.startswith("*"):
            name = m.group(1).upper()
            params.setdefault(name, no)
            nval = len([v for v in m.group(3).split(",")])
            if m.group(2) and nval > 10:
                issues.append((no, f"{name}(...) = {nval} values (max 10 per line)"))
        if cmd in ("*DIM", "*SET", "*GET", "*VGET") and len(f) > 1:
            params.setdefault(re.sub(r"\(.*", "", f[1]).upper(), no)
        if cmd == "*DO" and len(f) > 1:
            params.setdefault(f[1].upper(), no)
        if cmd == "CM" and len(f) > 1:
            comps.setdefault(f[1].upper(), no)
        if cmd == "MPDATA" and len(f) - 4 > 6:
            issues.append((no, f"MPDATA with {len(f) - 4} values (max 6)"))
        if cmd == "TBDATA" and len(f) - 2 > 6:
            issues.append((no, f"TBDATA with {len(f) - 2} values (max 6)"))
        if cmd == "*VWRITE" and no < len(lines):
            fmt = lines[no].strip()
            if len(fmt) > 80:
                issues.append((no + 1, f"*VWRITE format is {len(fmt)} characters (keep <= 80)"))
    for n in params:
        ALL_PARAMS.setdefault(n, path)
    for n in comps:
        ALL_COMPS.setdefault(n, path)
    for kind, no in stack:
        issues.append((no, f"{kind} never closed"))
    for name, no in params.items():
        if name.startswith("_") or name.startswith("%"):
            continue
        if len(name) > 32:
            issues.append((no, f"parameter {name} longer than 32 characters"))
        if name in RESERVED:
            issues.append((no, f"parameter {name} clashes with an APDL function/label"))
        if name in comps:
            issues.append((no, f"parameter {name} has the same name as component (CM line {comps[name]})"))
    return issues


if __name__ == "__main__":
    total = 0
    for path in sys.argv[1:]:
        found = lint(path)
        total += len(found)
        print(f"{path}: {'OK' if not found else str(len(found)) + ' issue(s)'}")
        for no, msg in sorted(found):
            print(f"   line {no}: {msg}")
    # a RESUME brings back every parameter of the earlier macros, so a
    # component in one macro must not reuse a parameter name of another
    for name, where in ALL_COMPS.items():
        if name in ALL_PARAMS and ALL_PARAMS[name] != where:
            total += 1
            print(f"cross-file: component {name} ({where}) = parameter in {ALL_PARAMS[name]}")
    sys.exit(1 if total else 0)
