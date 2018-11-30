#One-time use of Tableize_wrapper and File_prep_for_R to compare Paralogue annotation with SIFT and REVEL scores. Not intended to be part of the pipeline

import os, sys, subprocess, re, shlex, subprocess, codecs

input_file1 = sys.argv[1] #Tableized file
input_file2 = sys.argv[2] #paralogs file

os.system("python /data/Share/nick/Paralog_Anno/loftee/src/tableize_vcf.py --vcf " + input_file2 + " --out " + input_file2 + "_tableized_org --do_not_minrep --include_id --vep_info SYMBOL,Protein_position,Amino_acids,Codons,Paralogue_Vars,SIFT_score,SIFT_score --split_by_transcript --canonical_only")
