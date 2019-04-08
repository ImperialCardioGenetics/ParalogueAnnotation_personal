import  os, sys, subprocess, re

#SCRIPT FOR CSV TO VCF
# input_file = sys.argv[1]	#path of HCM csv file e.g. /media/nick/Data/Users/N/Documents/PhD/Paralogues/data_files/case_controls/HCM_missense_LMM_OMGL_chr.csv
# dir1 = input_file.rsplit("/", 1)[0]
# input_file2 = sys.argv[2]	#path of ExAC csv file e.g. /media/nick/Data/Users/N/Documents/PhD/Paralogues/data_files/case_controls/ExAC_missense_sarcomeric_chr.csv

# out_file = open(input_file.rsplit(".", 1)[0]+".vcf", "w")
# out_file.write("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n")
# with open(input_file) as f:
# 	for line in f:
# 		if line[0].isdigit():
# 			line = line.rstrip().split(",")
# 			ref = line[2].split(">")[0][-1]
# 			alt = line[2].split(">")[1]
# 			out_file.write(
# 				str(line[8])+"\t"+
# 				str(line[9])+"\t"+
# 				str(line[0])+"\t"+
# 				str(ref)+"\t"+
# 				str(alt)+"\t.\t.\t"+
# 				# "gene="+str(line[4])+",HGVSc="+str(line[5])+",HGVSp="+str(line[6])+",pathogenic="+str(line[7])+"\n"
# 				".\n"
# 				)			
# out_file.close()

# out_file2 = open(input_file2.rsplit(".", 1)[0]+".vcf", "w")
# out_file2.write("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n")
# with open(input_file2) as f:
# 	for line in f:
# 		if line[0].isdigit():
# 			line = line.rstrip().split(",")
# 			ref = line[2].split(">")[0][-1]
# 			alt = line[2].split(">")[1]
# 			out_file2.write(
# 				str(line[7])+"\t"+
# 				str(line[8])+"\t"+
# 				str(line[0])+"\t"+
# 				str(ref)+"\t"+
# 				str(alt)+"\t.\t.\t"+
# 				# "gene="+str(line[4])+",HGVSc="+str(line[5])+",HGVSp="+str(line[6])+",pathogenic="+str(line[7])+"\n"
# 				".\n"
# 				)			
# out_file2.close()

input_file = sys.argv[1]	#path of ICC MUTATION csv file to convert to vcf e.g. /media/nick/Data/PhD/Paralogues/ParalogueAnnotation_personal/data/LQTS/Mayo_LQTS_variants_for_vcf_convertion.csv
dir1 = input_file.rsplit("/", 1)[0]

out_file = open(input_file.rsplit(".", 1)[0]+".vcf", "w")
out_file.write("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n")

with open(input_file) as f:
	for line in f:
		print(line)
		if line[0].isdigit():
			line = line.rstrip().split(",")
			chrom = line[0]
			pos = line[1]
			ID = line[2]
			ref = line[3].split(">")[0][-1]
			alt = line[3].split(">")[1]
			out_file.write(
				str(chrom)+"\t"+
				str(pos)+"\t"+
				str(ID)+"\t"+
				str(ref)+"\t"+
				str(alt)+"\t.\t.\t"+
				# "gene="+str(line[4])+",HGVSc="+str(line[5])+",HGVSp="+str(line[6])+",pathogenic="+str(line[7])+"\n"
				".\n"
				)			
out_file.close()

input_file = sys.argv[2]	#path of gnomad csv file to convert to vcf e.g. /media/nick/Data/PhD/Paralogues/ParalogueAnnotation_personal/data/LQTS/gnomAD_v2.1.1_ENSG00000180509_2019_03_31_15_09_18_KCNE1.csv
dir1 = input_file.rsplit("/", 1)[0]

out_file = open(input_file.rsplit(".", 1)[0]+".vcf", "w")
out_file.write("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n")

with open(input_file) as f:
	for line in f:
		print(line)
		if line[0].isdigit():
			line = line.rstrip().split(",")
			chrom = line[0]
			pos = line[1]
			ID = line[2]
			ref = line[3].split(">")[0][-1]
			alt = line[3].split(">")[1]
			out_file.write(
				str(chrom)+"\t"+
				str(pos)+"\t"+
				str(ID)+"\t"+
				str(ref)+"\t"+
				str(alt)+"\t.\t.\t"+
				# "gene="+str(line[4])+",HGVSc="+str(line[5])+",HGVSp="+str(line[6])+",pathogenic="+str(line[7])+"\n"
				".\n"
				)			
out_file.close()

