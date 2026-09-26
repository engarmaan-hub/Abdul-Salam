"""
Rebuild the geometry of ansys/01_thermal_transient.mac in Python and check it
without ANSYS:

  * replays the K, L (incl. *DO loops), FLST/FITEM/AL commands
    (ANSYS does not create a second line between the same two keypoints),
  * checks that every area is a closed loop of 4 lines,
  * prints the lines that carry boundary conditions and contact pairs,
  * finds every exterior line and whether it has a thermal BC,
  * checks mapped-mesh compatibility (MSHKEY,1 needs equal divisions on
    opposite sides) and estimates the element count.

Usage:  python3 tools/check_geometry.py [path/to/01_thermal_transient.mac]
"""
import math
import re
import sys

MAC = sys.argv[1] if len(sys.argv) > 1 else "ansys/01_thermal_transient.mac"


def strip(line):
    return line.split("!")[0].strip()


def evaluate(expr, par):
    expr = expr.strip()
    if expr == "":
        return 0.0
    for k, v in par.items():
        expr = re.sub(rf"\b{k}\b", repr(v), expr, flags=re.I)
    return float(eval(expr))


def parse(path):
    """Returns par, kp, lines, areas, lesize.  Also sets parse.ndiv (lines
    with an explicit, hard number of divisions) and parse.free (areas meshed
    with MSHKEY,0)."""
    par, kp, lines, areas, lesize = {}, {}, {}, [], {}
    ndiv, free, mshkey, asel = {}, set(), 1, []
    src = [strip(l) for l in open(path)]
    flst, i = [], 0

    def add_line(a, b):
        for n, (p, q) in lines.items():
            if {p, q} == {a, b}:
                return n                     # ANSYS reuses the existing line
        n = len(lines) + 1
        lines[n] = (a, b)
        return n

    while i < len(src):
        s = src[i]
        up = s.upper()
        m = re.match(r"^([A-Z_][A-Z0-9_]*)\s*=\s*([-+.\deE*/()A-Z_ ]+)$", s, re.I)
        if m and not up.startswith("*"):
            try:
                par[m.group(1).upper()] = evaluate(m.group(2), par)
            except Exception:
                pass
        f = [x.strip() for x in s.split(",")]
        cmd = f[0].upper()
        if cmd == "K":
            n = int(f[1])
            x = evaluate(f[2] if len(f) > 2 else "", par)
            y = evaluate(f[3] if len(f) > 3 else "", par)
            kp[n] = (x, y)
        elif cmd == "L":
            add_line(int(evaluate(f[1], par)), int(evaluate(f[2], par)))
        elif cmd == "*DO":
            var = f[1].upper()
            try:
                a, b = int(evaluate(f[2], par)), int(evaluate(f[3], par))
            except Exception:              # loops over run-time values: skip
                a, b = 1, 0
            body = []
            i += 1
            while src[i].upper() != "*ENDDO":
                body.append(src[i])
                i += 1
            for v in range(a, b + 1):
                par[var] = v
                for bl in body:
                    g = [x.strip() for x in bl.split(",")]
                    if g[0].upper() == "L":
                        add_line(int(evaluate(g[1], par)), int(evaluate(g[2], par)))
        elif cmd == "FLST":
            flst = []
        elif cmd == "FITEM":
            flst.append(int(f[2]))
        elif cmd == "AL":
            areas.append(list(flst))
        elif cmd == "MSHKEY":
            mshkey = int(f[1])
        elif cmd == "ASEL" and len(f) > 5 and f[1].upper() in ("S", "A", "") and f[4].isdigit():
            rng = list(range(int(f[4]), int(f[5]) + 1)) if len(f) > 5 and f[5] else [int(f[4])]
            asel = rng if f[1].upper() == "S" else asel + rng
        elif cmd == "AMESH" and len(f) > 1 and f[1].upper() == "ALL" and mshkey == 0:
            free |= set(asel)
        elif cmd == "LESIZE" and f[1].isdigit():
            ndiv[int(f[1])] = int(f[4])                 # LESIZE,line,,,NDIV,...
        elif cmd == "LESIZE":
            size = float(f[2])
            items = []
            for k, v in enumerate(flst):
                if v < 0:
                    items += list(range(flst[k - 1] + 1, -v + 1))
                else:
                    items.append(v)
            for n in items:
                lesize.setdefault(n, size)
        i += 1
    parse.ndiv, parse.free = ndiv, free
    return par, kp, lines, areas, lesize


def main():
    par, kp, lines, areas, lesize = parse(MAC)
    L = lambda n: (kp[lines[n][0]], kp[lines[n][1]])
    length = lambda n: math.dist(*L(n))
    fmt = lambda p: f"({p[0]:.3f},{p[1]:.3f})"

    print(f"parameters: XL={par.get('XL')}  YB={par.get('YB')}")
    print(f"{len(kp)} keypoints, {len(lines)} lines, {len(areas)} areas")
    ok = True

    # closed loops
    for a, ls in enumerate(areas, 1):
        ends = [lines[n] for n in ls]
        deg = {}
        for p, q in ends:
            deg[p] = deg.get(p, 0) + 1
            deg[q] = deg.get(q, 0) + 1
        if len(ls) != 4 or any(v != 2 for v in deg.values()):
            ok = False
            print(f"  AREA {a} NOT CLOSED: lines {ls}")
    print("area loops: all closed" if ok else "area loops: PROBLEM")

    # BC and contact lines
    groups = {
        "natural ground surface  D,TEMP,%dibiao%": [27, 28],
        "exposed concrete        D,TEMP,%qiaotai%": [1, 2, 11, 12],
        "embankment top          D,TEMP,%lumian% + traffic": [55, 56],
        "contact behind  target (concrete)": [3, 4, 5, 6, 7],
        "contact behind  contact (soil)": [43, 44, 45, 46, 47],
        "contact front   target (concrete)": [13, 22, 23],
        "contact front   contact (soil)": [52, 53, 54],
        "contact footing target (concrete)": [8, 24, 25, 26],
        "contact footing contact (soil)": [48, 49, 50, 51],
    }
    print("\nlines used for BCs and contact")
    for name, ls in groups.items():
        segs = ", ".join(f"{n}:{fmt(L(n)[0])}-{fmt(L(n)[1])}" for n in ls)
        print(f"  {name}\n     {segs}")

    # contact pairs coincide?
    def seg_set(ls):
        return {tuple(sorted((tuple(round(c, 4) for c in L(n)[0]),
                              tuple(round(c, 4) for c in L(n)[1])))) for n in ls}
    for tgt, con in (([3, 4, 5, 6, 7], [43, 44, 45, 46, 47]), ([13, 22, 23], [52, 53, 54]),
                     ([8, 24, 25, 26], [48, 49, 50, 51])):
        same = seg_set(tgt) == seg_set(con)
        print(f"  pair {tgt} / {con}: {'coincident' if same else 'NOT coincident'}")

    # exterior lines
    use = {}
    for ls in areas:
        for n in ls:
            use[n] = use.get(n, 0) + 1
    geo = {n: seg_set([n]).pop() for n in lines}
    bc = set(sum((groups[k] for k in list(groups)[:3]), []))
    print("\nexterior lines (one area, no coincident twin)")
    for n in sorted(lines):
        twins = [m for m in lines if m != n and geo[m] == geo[n]]
        if use.get(n, 0) == 1 and not twins:
            (x1, y1), (x2, y2) = L(n)
            where = ("left end" if abs(x1 - par["XL"]) < 1e-6 and abs(x2 - par["XL"]) < 1e-6 else
                     "right end" if abs(x1 - 15) < 1e-6 and abs(x2 - 15) < 1e-6 else
                     "bottom" if abs(y1 - par["YB"]) < 1e-6 and abs(y2 - par["YB"]) < 1e-6 else "")
            tag = "BC" if n in bc else ("adiabatic side" if where in ("left end", "right end")
                                        else "HFLUX" if where == "bottom" else "NO BC !")
            print(f"  {n:3d} {fmt(L(n)[0])}-{fmt(L(n)[1])}  {tag}")

    # mapped mesh: divisions from LESIZE, then carried across opposite
    # sides of each area (as ANSYS does for lines without LESIZE)
    print("\nmapped mesh (divisions on opposite sides)")

    def pairs(ls):
        l0 = ls[0]
        opp0 = [n for n in ls[1:] if not set(lines[n]) & set(lines[l0])][0]
        rest = [n for n in ls if n not in (l0, opp0)]
        return (l0, opp0), tuple(rest)

    for rule, fn in (("ceil", math.ceil), ("nint", lambda v: max(1, round(v)))):
        div = {n: fn(length(n) / s - 1e-9) for n, s in lesize.items()}
        div.update(parse.ndiv)
        mapped = [ls for a, ls in enumerate(areas, 1) if a not in parse.free]
        changed = True
        while changed:
            changed = False
            for ls in mapped:
                for p, q in pairs(ls):
                    if p in div and q not in div:
                        div[q] = div[p]; changed = True
                    elif q in div and p not in div:
                        div[p] = div[q]; changed = True
        bad = [(a, p, div.get(p), q, div.get(q)) for a, ls in enumerate(areas, 1)
               if a not in parse.free for p, q in pairs(ls) if div.get(p) != div.get(q)]
        missing = sorted(n for n in lines if n not in div)
        nel = sum(div[pairs(ls)[0][0]] * div[pairs(ls)[1][0]] for ls in areas
                  if pairs(ls)[0][0] in div and pairs(ls)[1][0] in div)
        print(f"  LESIZE rounding = {rule}: "
              f"{'all opposite sides match' if not bad else 'MISMATCH (area, line, div, line, div): ' + str(bad)}"
              f"; no size: {missing or 'none'}; ~{nel} elements, ~{nel + len(areas) * 40} nodes")
    print(f"  free-meshed areas: {sorted(parse.free) or 'none'}; hard divisions: {parse.ndiv or 'none'}")
    # interface conformity: concrete line vs soil line, divisions must match
    iface = ((8, 48), (25, 49), (26, 50), (24, 51), (13, 54), (22, 52), (23, 53),
             (3, 43), (4, 44), (5, 45), (6, 46), (7, 47))
    div = {n: math.ceil(length(n) / s - 1e-9) for n, s in lesize.items()}
    div.update(parse.ndiv)
    bad_if = [(c, div.get(c), q, div.get(q)) for c, q in iface if div.get(c) != div.get(q)]
    print("  concrete/soil interface divisions: " + ("all match (nodes coincide)" if not bad_if
          else "MISMATCH (concrete line, div, soil line, div): " + str(bad_if)))
    print("  (soft LESIZE divisions, KYNDIV = 1, are only changed by the mesher when")
    print("   opposite sides of a mapped area disagree; with none left, the mesh")
    print("   ANSYS builds is the one checked here)")

if __name__ == "__main__":
    main()
