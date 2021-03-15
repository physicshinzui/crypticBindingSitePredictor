#!/bin/bash
set -Ceu

read -p "Type 0 or 1 : " num

if [ ${num} -eq 0 ]; then 
    crypticSitePredictor.py -r samples/sample1/ref.pdb -i samples/sample1/traj_aligned_100.xtc -o test_file
    echo 'diff...' 
    diff sasa_test_file.csv  samples/sample1/sasa_base.csv
    diff rsasa_test_file.csv samples/sample1/rsasa_base.csv
    
    crypticSitePredictor.py -r samples/sample2/ref.pdb -i samples/sample2/traj_aligned_100.xtc -o test_file
    echo 'diff...' 
    diff sasa_test_file.csv samples/sample2/sasa_base.csv
    diff rsasa_test_file.csv samples/sample2/rsasa_base.csv
 


elif [ ${num} -eq 1 ]; then 
    crypticSitePredictor.py -r samples/sample3/protein-h_renum.pdb -i samples/sample3/all_skip100.xtc -o test_bclxl_file
    echo 'diff...' 
    diff sasa_test_file.csv samples/sample3/sasa_base.csv
    diff rsasa_test_file.csv samples/sample3/rsasa_base.csv

elif [ ${num} -eq 2 ]; then 
    traj=samples/sample4/all_skip100.xtc
    ref=samples/sample4/ref_protein-h.pdb
    crypticSitePredictor.py -r ${ref} -i ${traj} -o test_integrin_file
fi
