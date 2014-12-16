File Formats
============

The formats of the .data and .cov files are documented at:

  http://darkmatter.ps.uci.edu/wiki/DeepZot/Baofit/Format

The .scan files contain three columns: alpha-perp, alpha-parallel, chisq


How to reproduce published BOSS results
=======================================

The following instructions assume that you are in a subdirectory of the baofit package and have installed the baofit program so that it is in your search path:

# Busca 2013
baofit -i ../config/BOSSDR9LyaFXi.ini

# Slosar 2013
baofit -i ../config/BOSSDR9LyaF.ini

# Font-Ribera 2013
baofit -i ../config/BOSSDR11QSOLyaF.ini

# Delubac 2014
baofit -i ../config/BOSSDR11LyaF_k.ini

All of the resulting output files will be created in your working directory and prefixed with XXX_ (with XXX = BOSSDR9LyaFXi, etc).


Chi-Square Scans
================

To reproduce the .scan files, add the --parameter-scan option to produce XXX_scan.dat. This file starts with a three-line header that describes the baseline fit, followed by one line per scan point. The scan point lines give each parameter value (including fixed parameters) in the same order that they are normally displayed, followed by the chisq value. The anisotropic BAO scale parameters, alpha-parallel and alpha-perp, are columns 6 and 7 (auto-correlation; r-space model), 8 and 9 (auto-correlation; k-space model) or 9 and 10 (cross-correlation; r-space model) counting from 0.

You can use the parsescan.py python script to extract the 3 columns that appear in the .scan files (note that we swap the order of the alphas here):

# Busca 2013
./parsescan.py BOSSDR9LyaFXi_scan.dat BOSSDR9LyaFXi.scan 7 6 

# Slosar 2013
./parsescan.py BOSSDR9LyaF_scan.dat BOSSDR9LyaF.scan 7 6

# Font-Fibera 2013
./parsescan.py BOSSDR11QSOLyaF_scan.dat BOSSDR11QSOLyaF.scan 10 9

# Delubac 2014
./parsescan.py BOSSDR11LyaF_k_scan.dat BOSSDR11LyaF_k.scan 9 8

There is a mathematica package for more sophisticated processing of baofit outputs:

  https://github.com/deepzot/mathpkg/blob/master/BaoFitTools.m

For example, the chisq scan plots above can be reproduced using, e.g.

Needs["DeepZot`BaoFitTools`"]
loadFitAnalysis[BOSSDR9LyaFXi, "BOSSDR9LyaFXi_scan"]
fitContoursPlot[{fitAnalysisSamples[BOSSDR9LyaFXi, parameters -> {8, 7}]}, 
  FrameLabel -> {"a(perp)","a(parallel)"}, styles -> {Thick}, Axes -> None]
