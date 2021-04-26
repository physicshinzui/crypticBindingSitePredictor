# Cryptic Binding Site Precictor
Author: Shinji Iida

# Requirements
MDAnalysis >= 1.0.0 (https://www.mdanalysis.org/)
mdtraj     >= 1.9.4 (https://www.mdtraj.org/1.9.4/index.html#)
pandas     >= 1.1.3
python     >= 3.6

# Installation
`export PATH=${PATH}:path_to_/crypticSitePredictor/src`

# Usage
```
$ crypticSitePredictor.py -h
usage: crypticSitePredictor.py [-h] -r REF -i TRJ [-o OUT] [-a ALPHA] [-b BETA] [-th THRESHOLD] [-tr TRJ_RANGE [TRJ_RANGE ...]]

This is Cryptic-site predictor

optional arguments:
  -h, --help            show this help message and exit
  -r REF, --ref REF
  -i TRJ, --trj TRJ
  -o OUT, --out OUT     output suffix
  -a ALPHA, --alpha ALPHA
                        Upper bound of Delta F
  -b BETA, --beta BETA  Lower bound of sigma
  -th THRESHOLD, --threshold THRESHOLD
                        threshold
  -tr TRJ_RANGE [TRJ_RANGE ...], --trj_range TRJ_RANGE [TRJ_RANGE ...]
                        traj range
```

# Notes
- More than 5 independent canonical MD simulations more than 500 ns are preferable for the prediction. 
- Input trajectories are assumed to have no hydrogen atom. 
- Golden section spiral algorithm calcuates Solvent Accessible Surface Area per residue

# References
Please cite:
https://dx.doi.org/10.1021/acs.jpcb.0c04963

# Licence 

