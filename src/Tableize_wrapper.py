import os, sys, subprocess, re

#Tableize.py wrapper plus additional formatting that can't be done in tableize. Also adds in Para Z scores.

def Tableize_wrap(input_file):
	#input_file = sys.argv[1]	#path of paraloc file e.g. /data/Share/nick/Paralog_Anno/data_files/clinvar_20171029_onlyPathogenic.out_paraloc
	dir1 = input_file.rsplit("/", 1)[0]

	os.system("python /data/Share/nick/Paralog_Anno/loftee/src/tableize_vcf.py --vcf " + input_file + " --out " + input_file + "_tableized_org --do_not_minrep --include_id --vep_info SYMBOL,Protein_position,Amino_acids,Codons,Paralogue_Vars --split_by_transcript --canonical_only")

	out_file = open(input_file + "_tableized", "w")
	with open(input_file + "_tableized_org") as f:
		for line in f:
			if line.startswith("CHROM"):
				out_file.write(line.rstrip()+"\tPara_Z_score\n")
			else:
				line = line.split()
				AAs = line[8].split(",")
				Codons = line[9].split(",")
				Gene = str(line[6])+".txt"
				for column in range(0,len(AAs)):
					if (
						line[7] != "NA" and 
						not "-" in line[7]
						):
						if Gene in os.listdir("/data/Share/nick/Paralog_Anno/data_files/para_zscores/genes"):
							with open("/data/Share/nick/Paralog_Anno/data_files/para_zscores/genes/"+Gene) as para_z_file:
								for i, para_z_line in enumerate(para_z_file):
									if i == int(line[7])-1:
										if AAs[column].split("/")[0] == para_z_line.split()[1]:
											para_z_score = para_z_line.split()[2]
										else:
											para_z_score = "NA"
									elif i > int(line[7])-1:
										break
						else:
							para_z_score = "NA"
					else:
						para_z_score = "NA"
					out_file.write(line[0] + "\t" + line[1] +"\t" +line[2] +"\t" + line[3] +"\t" +line[4] +"\t" +line[5] +"\t" + line[6] + "\t" + line[7] + "\t" + AAs[column] + "\t" + Codons[column] + "\t" + line[10] + "\t" + str(para_z_score) + "\n")
		out_file.close()

	print("Tableize_wrapper Done!")
