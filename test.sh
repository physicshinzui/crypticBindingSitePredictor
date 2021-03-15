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
 
    crypticSitePredictor.py -r samples/sample3/protein-h_renum.pdb -i samples/sample3/all_skip100.xtc -o test_file -a 5.0
    echo 'diff...' 
    diff sasa_test_file.csv samples/sample3/sasa_base.csv
    diff rsasa_test_file.csv samples/sample3/rsasa_base.csv

elif [ ${num} -eq 1 ]; then 
    traj=/Volumes/siida-SSD/cryptic_project/bcl-xl/processed_xtc/protein-h_frm200ns/npt_prod_2.xtc
    ref=/Volumes/siida-SSD/cryptic_project/bcl-xl/ref/protein-h_renum.pdb
    crypticSitePredictor.py -r ${ref} -i ${traj} -o test2_file

elif [ ${num} -eq 2 ]; then 
    traj=/Volumes/siida-SSD/cryptic_project/integrin/processed_xtc/protein-h_frm200ns/for_test/all_skip100.xtc
    ref=/Volumes/siida-SSD/cryptic_project/integrin/processed_xtc/protein-h_frm200ns/for_test/ref_protein-h.pdb
    crypticSitePredictor.py -r ${ref} -i ${traj} -o test_integrin_file
fi
