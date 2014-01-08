How to reproduce published BOSS results
=======================================

The following instructions assume that you are in a subdirectory of the baofit package and have installed the baofit program so that it is in your search path:

# Busca 2013
baofit -i ../config/BOSSDR9LyaFXi.ini

# Slosar 2013
baofit -i ../config/BOSSDR9LyaF.ini

# Font-Ribera 2013
baofit -i ../config/BOSSDR11QSOLyaF.ini

All of the resulting output files will be created in your working directory and prefixed with XXX_ (with XXX = BOSSDR9LyaFXi, etc).

File Formats
============

The formats of the XXX.data and XXX.cov files are documented at:

  http://darkmatter.ps.uci.edu/wiki/DeepZot/Baofit/Format

The XXX_scan.dat files contain three columns (alpha-parallel, alpha-perp, chisq) and can be reproduced by adding the --parameter-scan option to any of command-lines given above.
