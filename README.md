# configtools
Various tools to manipulate configurations of atoms/molecules for input to simulation codes

tip4p
-----

Usage : `tip4p <xmolfile> <style>`

Reads an xmol file containing (1) a list of cell vectors on line 2 and (2) a list of coordinates
for oxygen and hydrogen atoms in any order. Adjusts and outputs these to be consistent with the tip4p
rigid water model with style being either tip4p, tip4p/Ice or tip4p/2005. Output includes the coordinates
of the massless charge site in the H-O-H bisector.

xmol2config
-----------

Usage: `xmol2config <xmolfile>`

Reads specified xmol file and prints equivalent DL_POLY CONFIG file to stdout. Expects to find the matrix of cell vectors on the comment line of the input xmol file, listed as 9 components.                     

