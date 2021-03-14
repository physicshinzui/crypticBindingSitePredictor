#!/bin/bash
set -Ceu



crypticSitePredictor.py -r samples/sample1/ref.pdb -i samples/sample1/traj_aligned_100.xtc -o test_file
echo 'diff...' 
diff sasa_test_file.csv  samples/sample1/sasa_base.csv
diff rsasa_test_file.csv samples/sample1/rsasa_base.csv

crypticSitePredictor.py -r samples/sample2/ref.pdb -i samples/sample2/traj_aligned_100.xtc -o test_file
echo 'diff...' 
diff sasa_test_file.csv samples/sample2/sasa_base.csv
diff rsasa_test_file.csv samples/sample2/rsasa_base.csv
 