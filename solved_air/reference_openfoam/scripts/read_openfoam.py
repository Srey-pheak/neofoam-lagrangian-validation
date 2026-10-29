#!/usr/bin/env python3
"""Convert the lagrangian output of an OpenFOAM run to the NeoFOAM CSV layout.

Reads <time>/lagrangian/kinematicCloud/{positions,U,d,origId} of every time folder
(in v2406 'positions' holds Cartesian "(x y z) celli" entries) and writes
    t,id,x,y,z,ux,uy,uz,d          id = origId (injection order)

    python3 scripts/read_openfoam.py <run folder> <output csv>
"""
import re
import sys
from pathlib import Path

NUM = r"[-+]?(?:\d+\.?\d*|\.\d+)(?:[eE][-+]?\d+)?"
VEC = re.compile(r"\(\s*(" + NUM + r")\s+(" + NUM + r")\s+(" + NUM + r")\s*\)")


def body(path):
    """file text without the FoamFile header and comments"""
    text = path.read_text()
    i = text.find("FoamFile")
    if i >= 0:
        text = text[text.find("}", i) + 1:]
    return re.sub(r"//.*", "", text)


def read_vectors(path):
    text = body(path)
    m = re.match(r"\s*(\d+)\s*\{\s*(\(.*?\))\s*\}", text, re.S)   # uniform list  N{(x y z)}
    if m:
        v = tuple(float(c) for c in VEC.search(m.group(2)).groups())
        return [v] * int(m.group(1))
    return [tuple(float(c) for c in m.groups()) for m in VEC.finditer(text)]


def read_scalars(path):
    text = body(path)
    m = re.match(r"\s*(\d+)\s*\{\s*(" + NUM + r")\s*\}", text)       # uniform list  N{v}
    if m:
        return [float(m.group(2))] * int(m.group(1))
    m = re.match(r"\s*(\d+)\s*\((.*)\)", text, re.S)                  # N( v0 v1 ... )
    return [float(v) for v in m.group(2).split()]


def read_U(path):
    """cell values of a volVectorField file (internalField, uniform or nonuniform) -> list of (x y z)"""
    text = body(path)
    field = text[text.find("internalField"):text.find("boundaryField")]
    if "nonuniform" not in field:
        return None                                                    # uniform: caller knows the mesh
    return [tuple(float(c) for c in m.groups()) for m in VEC.finditer(field)]


def main():
    case, out = Path(sys.argv[1]), Path(sys.argv[2])
    rows = []
    for d in case.iterdir():
        try:
            t = float(d.name)
        except ValueError:
            continue
        c = d / "lagrangian" / "kinematicCloud"
        if not (c / "positions").exists():
            continue
        xs, us = read_vectors(c / "positions"), read_vectors(c / "U")
        ds, ids = read_scalars(c / "d"), read_scalars(c / "origId")
        for i, x, u, dd in sorted(zip(ids, xs, us, ds)):
            rows.append((t, int(i), *x, *u, dd))
    rows.sort()
    if not rows:
        sys.exit(f"no lagrangian output in {case} (not run yet?)")
    out.parent.mkdir(parents=True, exist_ok=True)
    with open(out, "w") as f:
        f.write("t,id,x,y,z,ux,uy,uz,d\n")
        for r in rows:
            f.write(",".join(repr(v) for v in r) + "\n")
    print(f"{len(rows)} rows ({len({r[0] for r in rows})} times) -> {out}")


if __name__ == "__main__":
    main()
