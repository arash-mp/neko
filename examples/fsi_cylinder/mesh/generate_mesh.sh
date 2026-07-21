#!/bin/bash

# Set the paths to your executables here

GMSH_PATH="/Gmsh-Path/bin/gmsh"
GMSH2NEK_PATH="/Neko-Path/contrib/gmsh2nek/gmsh2nek"
REA2NBIN_PATH="/Neko-Path/bin/rea2nbin"

"$GMSH_PATH" cylinder.geo -

"$GMSH2NEK_PATH" <<EOF
3
3D_ext_cyl
1
4 3
0.0 0.0 1.0
3D_ext_cyl
EOF

"$REA2NBIN_PATH" 3D_ext_cyl.re2 fsi_cylinder.nmsh
cp fsi_cylinder.nmsh ../

echo "Process completed successfully."

