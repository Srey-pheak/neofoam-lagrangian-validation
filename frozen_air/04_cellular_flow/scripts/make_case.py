#!/usr/bin/env python3
"""Test 4 - writes the case files that depend on the parameters (scripts/params.py):

  system/blockMeshDict                     2L x 2L x dz box, nx x ny x 1 cells
  system/controlDict                       dt, endTime, writeInterval
  system/setExprFieldsDict                 the cellular-flow formula (setExprFields writes it into 0/U)
  0.orig/U                                 start field (zero; overwritten by setExprFields)
  constant/kinematicCloudProperties        stokesDrag, no gravity, RK2, findCell, one injector per St
  constant/kinematicCloudPositions.St<St>  the start positions of each size

Run once from the test folder after changing params.py:   python3 scripts/make_case.py
"""
import math
import os
import sys

sys.path.insert(0, os.path.dirname(__file__))
from params import *  # noqa: E402,F403

HEADER = """/*--------------------------------*- C++ -*----------------------------------*\\
| =========                 |                                                 |
| \\\\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox           |
|  \\\\    /   O peration     | Version:  v2406                                 |
|   \\\\  /    A nd           | Website:  www.openfoam.com                      |
|    \\\\/     M anipulation  |                                                 |
\\*---------------------------------------------------------------------------*/
FoamFile
{{
    version     2.0;
    format      ascii;
    class       {cls};
    object      {obj};
}}
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
// written by scripts/make_case.py from scripts/params.py: change params.py, not this file

"""
END = "\n// ************************************************************************* //\n"


def write(path, cls, obj, body):
    with open(path, "w") as f:
        f.write(HEADER.format(cls=cls, obj=obj) + body + END)
    print("written", path)


def st_name(St):
    return f"St{St:g}"


# ─── mesh ───
write("system/blockMeshDict", "dictionary", "blockMeshDict", f"""scale   1;

Lx      {2 * L:.17g};
Lz      {dz:.17g};

vertices
(
    (0 0 0) ($Lx 0 0) ($Lx $Lx 0) (0 $Lx 0)
    (0 0 $Lz) ($Lx 0 $Lz) ($Lx $Lx $Lz) (0 $Lx $Lz)
);

blocks
(
    hex (0 1 2 3 4 5 6 7) ({nx} {ny} 1) simpleGrading (1 1 1)
);

boundary
(
    sides
    {{
        type patch;
        faces
        (
            (0 3 2 1) (4 5 6 7) (0 4 7 3) (1 2 6 5) (0 1 5 4) (3 7 6 2)
        );
    }}
);
""")

# ─── time ───
write("system/controlDict", "dictionary", "controlDict", f"""application     neoIcoUncoupledKinematicParcelFoam;

startFrom       startTime;
startTime       0;
stopAt          endTime;
endTime         {end_time:g};
deltaT          {dt:g};

writeControl    timeStep;
writeInterval   {write_interval};             // particles.csv every {write_interval} steps
purgeWrite      0;
writeFormat     ascii;
writePrecision  16;
writeCompression off;
timeFormat      general;
timePrecision   10;
runTimeModifiable false;
""")

# ─── cellular flow: the formula in system/setExprFieldsDict, evaluated by setExprFields
#     at every cell centre (pos()) before the solver runs; 0.orig/U is only the start field ───
write("system/setExprFieldsDict", "dictionary", "setExprFieldsDict", f"""expressions
(
    U
    {{
        field       U;
        dimensions  [0 1 -1 0 0 0 0];

        variables
        (
            "L = {L:.17g}"
            "U0 = {U0:.17g}"
        );

        // cellular flow (Maxey 1987): 2 x 2 vortices of size L
        expression
        #{{
            vector
            (
                 U0*sin(pi()*pos().x()/L)*cos(pi()*pos().y()/L),
                -U0*cos(pi()*pos().x()/L)*sin(pi()*pos().y()/L),
                0
            )
        #}};
    }}
);
""")

os.makedirs("0.orig", exist_ok=True)
write("0.orig/U", "volVectorField", "U", """dimensions      [0 1 -1 0 0 0 0];

// start field: setExprFields (system/setExprFieldsDict) writes the cellular flow into 0/U
internalField   uniform (0 0 0);

boundaryField
{
    sides
    {
        type            zeroGradient;
    }
}
""")

# ─── particles ───
pos = start_positions()
injectors = []
for St in ST:
    d = diameter(St)
    name = st_name(St)
    write(f"constant/kinematicCloudPositions.{name}", "vectorField", f"kinematicCloudPositions.{name}",
          f"{len(pos)}\n(\n" + "\n".join(f"({x:.17g} {y:.17g} {z:.17g})" for x, y, z in pos) + "\n)\n")
    injectors.append(f"""        particle_{name}        // St = {St:g}: d = {d * 1e6:.2f} um, tau_p = {St * tau_f:g} s
        {{
            type            manualInjection;
            massTotal       0;
            parcelBasisType fixed;
            nParticle       1;
            SOI             0;
            positionsFile   "kinematicCloudPositions.{name}";
            U0              (0 0 0);
            sizeDistribution
            {{
                type        fixedValue;
                fixedValueDistribution
                {{
                    value   {d:.17g};
                }}
            }}
        }}""")

write("constant/kinematicCloudProperties", "dictionary", "kinematicCloudProperties", """solution
{
    active          true;
    coupled         false;      // one-way: air is frozen
    transient       yes;
    cellValueSourceCorrection off;
    maxCo           1.0;

    locate          findCell;   // OpenFOAM findCell, cell value of U

    interpolationSchemes
    {
        rho             cell;
        U               cell;
        mu              cell;
    }

    integrationSchemes
    {
        U               RK2;
    }
}

constantProperties
{
    rho0            2500;       // solid particle
    youngsModulus   1e9;        // required by colliding cloud, unused (no collisions)
    poissonsRatio   0.35;
}

subModels
{
    particleForces
    {
        stokesDrag;             // linear drag, as Maxey (1987); no gravity
    }

    injectionModels
    {
""" + "\n\n".join(injectors) + """
    }

    dispersionModel          none;
    patchInteractionModel    standardWallInteraction;
    standardWallInteractionCoeffs
    {
        type        rebound;            // unused: no walls
    }
    surfaceFilmModel         none;
    stochasticCollisionModel none;
    collisionModel           none;
}

cloudFunctions
{}
""")

print(f"\nSt -> d:  " + ",  ".join(f"St {St:g}: {diameter(St) * 1e6:.1f} um" for St in ST))
print(f"particles: {len(ST)} x {len(pos)};  max Re ~ {rho_f * U0 * diameter(max(ST)) / mu:.2f}")
