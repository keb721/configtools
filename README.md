# configtools
Various tools to manipulate configurations of atoms/molecules for input to simulation codes

tip4p
-----

Usage : `tip4p <xmolfile> <tip4p_style> <output_style>`

Reads an xmol file containing (1) a list of cell vectors on line 2 and (2) a list of coordinates
for oxygen and hydrogen atoms in any order. Adjusts and outputs these to be consistent with the tip4p
rigid water model with style being either tip4p, tip4p/Ice, tip4p/2005 or tip4p/long (recommended 
parameters for use with long range Coulombic solver). This can output for either dlpoly, gromacs, 
lammpsm (output includes the coordinates of the massless charge site in the H-O-H bisector), or 
lammps (output only includes oxygens and hydrogens). By default the gromacs topolgy file is written as 
``tip4p_gromacs.top''.

xmol2config
-----------

Usage: `xmol2config <xmolfile>`

Reads specified xmol file and prints equivalent DL_POLY CONFIG file to stdout. Expects to find the matrix of cell vectors on the comment line of the input xmol file, listed as 9 components.                     

