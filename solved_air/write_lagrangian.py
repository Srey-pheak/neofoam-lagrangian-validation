#!/usr/bin/env python3
"""Write the NeoFOAM particle CSV of a run as OpenFOAM lagrangian cloud files, for ParaView.

NeoFOAM writes the particles only to postProcessing/lagrangian/<cloud>/particles.csv. This script
writes them into the time folders the way OpenFOAM does,
    <time>/lagrangian/<cloud>/{positions, U, d, origId}
so ParaView's OpenFOAM reader (open the run's .foam file) shows the NeoFOAM particles exactly as
it shows OpenFOAM's. 'positions' is the "(x y z) celli" layout OpenFOAM v2406 writes; the cell is
the one whose centre is nearest (OpenFOAM postProcess -func writeCellCentres), which is the
containing cell on these hex meshes. <time>/uniform/lagrangian/<cloud>/cloudProperties says
"geometry positions", so OpenFOAM tools (foamToVTK, paraFoam) read 'positions' instead of the
barycentric 'coordinates' file OpenFOAM itself writes.

Needs OpenFOAM loaded (for writeCellCentres). Usage:
    python3 write_lagrangian.py <run folder> [cloud name, default kinematicCloud]
"""
import re
import shutil
import subprocess
import sys
from pathlib import Path

import numpy as np

NUM = r"[-+]?(?:\d+\.?\d*|\.\d+)(?:[eE][-+]?\d+)?"
VEC = re.compile(r"\(\s*(" + NUM + r")\s+(" + NUM + r")\s+(" + NUM + r")\s*\)")

HEADER = """FoamFile
{{
    version     2.0;
    format      ascii;
    arch        "LSB;label=32;scalar=64";
    class       {cls};
    location    "{loc}";
    object      {obj};
}}
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

"""


def cell_centres(run):
    """cell centres of the mesh, via OpenFOAM writeCellCentres (files removed afterwards)"""
    subprocess.run(["postProcess", "-func", "writeCellCentres", "-time", "0"], cwd=run, check=True,
                   stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
    text = (run / "0" / "C").read_text()
    field = text[text.find("internalField"):text.find("boundaryField")]
    C = np.array([[float(c) for c in m.groups()] for m in VEC.finditer(field)])
    for f in ("C", "Cx", "Cy", "Cz", "V"):
        (run / "0" / f).unlink(missing_ok=True)
    return C


def write(path, cls, obj, loc, body):
    path.write_text(HEADER.format(cls=cls, loc=loc, obj=obj) + body
                    + "\n\n// ************************************************************************* //\n")


def main():
    run = Path(sys.argv[1]).resolve()
    cloud = sys.argv[2] if len(sys.argv) > 2 else "kinematicCloud"
    csv = run / "postProcessing" / "lagrangian" / cloud / "particles.csv"
    if not csv.is_file():
        sys.exit(f"no particle CSV in {run}")
    D = np.genfromtxt(csv, delimiter=",", names=True)
    C = cell_centres(run)

    # the time folders the solver wrote (the CSV times carry round-off: match to the nearest folder)
    folders = {}
    for d in run.iterdir():
        try:
            folders[float(d.name)] = d
        except ValueError:
            pass
    ft = np.array(sorted(folders))

    n_written = 0
    for t in np.unique(D["t"]):
        k = np.argmin(np.abs(ft - t))
        if abs(ft[k] - t) > 1e-9 * max(1.0, abs(t)):
            continue                                      # no time folder for this output
        folder = folders[ft[k]]
        P = D[D["t"] == t]
        P = P[np.argsort(P["id"])]
        x = np.stack([P["x"], P["y"], P["z"]], axis=1)
        cells = np.argmin(((x[:, None, :] - C[None, :, :]) ** 2).sum(axis=2), axis=1)

        out = folder / "lagrangian" / cloud
        if out.exists():
            shutil.rmtree(out)
        out.mkdir(parents=True)
        loc = f"{folder.name}/lagrangian/{cloud}"
        n = len(P)
        write(out / "positions", "Cloud<basicKinematicCollidingParcel>", "positions", loc,
              f"{n}\n(\n" + "".join(f"({float(a)!r} {float(b)!r} {float(c)!r}) {int(ci)}\n" for (a, b, c), ci in zip(x, cells)) + ")")
        write(out / "U", "vectorField", "U", loc,
              f"{n}\n(\n" + "".join(f"({float(a)!r} {float(b)!r} {float(c)!r})\n" for a, b, c in zip(P["ux"], P["uy"], P["uz"])) + ")")
        write(out / "d", "scalarField", "d", loc, f"{n}\n(\n" + "".join(f"{float(v)!r}\n" for v in P["d"]) + ")")
        write(out / "origId", "labelField", "origId", loc,
              f"{n}\n(\n" + "".join(f"{int(v)}\n" for v in P["id"]) + ")")
        uni = folder / "uniform" / "lagrangian" / cloud
        uni.mkdir(parents=True, exist_ok=True)
        write(uni / "cloudProperties", "dictionary", "cloudProperties", f"{folder.name}/uniform/lagrangian/{cloud}",
              f"geometry        positions;\n\nprocessor0\n{{\n    particleCount   {n};\n}}")
        n_written += 1

    print(f"{n_written} time folders: {cloud} particles written as <time>/lagrangian/{cloud}/ in {run}")


if __name__ == "__main__":
    main()
