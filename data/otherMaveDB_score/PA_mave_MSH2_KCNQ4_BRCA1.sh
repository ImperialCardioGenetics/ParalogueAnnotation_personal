module load anaconda3/personal
source /home/nyl112/.bashrc

python3 ../../src/Paralogue_Annotation_masterscript_afterrun.py MSH2_mave_deleterious_b38.out_paraloc 38 paraloc noQC /work/nyl112/loftee/src/ /work/nyl112/data
python3 ../../src/Paralogue_Annotation_masterscript_afterrun.py MSH2_mave_benign_b38.out_paraloc 38 paraloc noQC /work/nyl112/loftee/src/ /work/nyl112/data
python3 ../../src/Paralogue_Annotation_masterscript_afterrun.py MSH2_mave_deleterious_b38.out_paraloc 38 paraloc para_con /work/nyl112/loftee/src/ /work/nyl112/data
python3 ../../src/Paralogue_Annotation_masterscript_afterrun.py MSH2_mave_benign_b38.out_paraloc 38 paraloc para_con /work/nyl112/loftee/src/ /work/nyl112/data
python3 ../../src/Paralogue_Annotation_masterscript_afterrun.py MSH2_mave_deleterious_b38.out_paraloc 38 paraloc all_con /work/nyl112/loftee/src/ /work/nyl112/data
python3 ../../src/Paralogue_Annotation_masterscript_afterrun.py MSH2_mave_benign_b38.out_paraloc 38 paraloc all_con /work/nyl112/loftee/src/ /work/nyl112/data

python3 ../../src/Paralogue_Annotation_masterscript_afterrun.py KCNQ4_mave_deleterious_b38.out_paraloc 38 paraloc noQC /work/nyl112/loftee/src/ /work/nyl112/data
python3 ../../src/Paralogue_Annotation_masterscript_afterrun.py KCNQ4_mave_benign_b38.out_paraloc 38 paraloc noQC /work/nyl112/loftee/src/ /work/nyl112/data
python3 ../../src/Paralogue_Annotation_masterscript_afterrun.py KCNQ4_mave_deleterious_b38.out_paraloc 38 paraloc para_con /work/nyl112/loftee/src/ /work/nyl112/data
python3 ../../src/Paralogue_Annotation_masterscript_afterrun.py KCNQ4_mave_benign_b38.out_paraloc 38 paraloc para_con /work/nyl112/loftee/src/ /work/nyl112/data
python3 ../../src/Paralogue_Annotation_masterscript_afterrun.py KCNQ4_mave_deleterious_b38.out_paraloc 38 paraloc all_con /work/nyl112/loftee/src/ /work/nyl112/data
python3 ../../src/Paralogue_Annotation_masterscript_afterrun.py KCNQ4_mave_benign_b38.out_paraloc 38 paraloc all_con /work/nyl112/loftee/src/ /work/nyl112/data

python3 ../../src/Paralogue_Annotation_masterscript_afterrun.py BRCA1_mave_patho_b38.out_paraloc 38 paraloc noQC /work/nyl112/loftee/src/ /work/nyl112/data
python3 ../../src/Paralogue_Annotation_masterscript_afterrun.py BRCA1_mave_benign_b38.out_paraloc 38 paraloc noQC /work/nyl112/loftee/src/ /work/nyl112/data
python3 ../../src/Paralogue_Annotation_masterscript_afterrun.py BRCA1_mave_patho_b38.out_paraloc 38 paraloc para_con /work/nyl112/loftee/src/ /work/nyl112/data
python3 ../../src/Paralogue_Annotation_masterscript_afterrun.py BRCA1_mave_benign_b38.out_paraloc 38 paraloc para_con /work/nyl112/loftee/src/ /work/nyl112/data
python3 ../../src/Paralogue_Annotation_masterscript_afterrun.py BRCA1_mave_patho_b38.out_paraloc 38 paraloc all_con /work/nyl112/loftee/src/ /work/nyl112/data
python3 ../../src/Paralogue_Annotation_masterscript_afterrun.py BRCA1_mave_benign_b38.out_paraloc 38 paraloc all_con /work/nyl112/loftee/src/ /work/nyl112/data
