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
		# print(line)
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

out_file2 = open(input_file.rsplit(".", 1)[0]+"_rare.vcf", "w")
out_file2.write("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n")

out_file3 = open(input_file.rsplit(".", 1)[0]+"_with_customIDs.csv", "w")

try:
	gene = sys.argv[3]
except IndexError:
	pass

with open(input_file) as f:
	ID_no = 1
	header_out_file3_check = 0
	for og_line in f:
		# print(line)
		if header_out_file3_check == 0:
			out_file3.write("Gene,"+og_line)
			header_out_file3_check = 1
		if og_line[0].isdigit():
			line = og_line.rstrip().split(",")
			csv_line = og_line.split(",",3)
			chrom = line[0]
			pos = line[1]
			# if line[2]: 
			# 	ID = line[2]
			# elif not line[2]:	#TEST THIS 
			# 	ID = "custom_" + str(ID_no)
			# 	ID_no += 1
			ID = "custom_" + str(ID_no)
			ID_no += 1
			ref = line[3]
			alt = line[4]
			out_file.write(
				str(chrom)+"\t"+
				str(pos)+"\t"+
				str(ID)+"\t"+
				str(ref)+"\t"+
				str(alt)+"\t.\t.\t"+
				# "gene="+str(line[4])+",HGVSc="+str(line[5])+",HGVSp="+str(line[6])+",pathogenic="+str(line[7])+"\n"
				".\n"
				)		
			allele_freq = float(line[15])
			if allele_freq <= 8.2e-06:
				out_file2.write(
					str(chrom)+"\t"+
					str(pos)+"\t"+
					str(ID)+"\t"+
					str(ref)+"\t"+
					str(alt)+"\t.\t.\t"+
					# "gene="+str(line[4])+",HGVSc="+str(line[5])+",HGVSp="+str(line[6])+",pathogenic="+str(line[7])+"\n"
					".\n"
					)
			out_file3.write(
				str(gene)+","+
				str(chrom)+","+
				str(pos)+","+
				str(ID)+","+
				str(csv_line[3])
				)

out_file.close()
out_file2.close()
out_file3.close()

'''
input_file = sys.argv[1]	#path of ICC MUTATION csv file to convert to vcf e.g. /media/nick/Data/PhD/Paralogues/ParalogueAnnotation_personal/data/LQTS/Mayo_LQTS_variants_for_vcf_convertion.csv
disease_cohort = input_file.rsplit(".", 1)[0]+".vcf"
input_file2 = sys.argv[2]	#path of gnomad csv file to convert to vcf e.g. /media/nick/Data/PhD/Paralogues/ParalogueAnnotation_personal/data/LQTS/gnomAD_v2.1.1_ENSG00000180509_2019_03_31_15_09_18_KCNE1.csv
gnomad_cohort = input_file2.rsplit(".", 1)[0]+".vcf"
gnomad_cohort_rare = input_file2.rsplit(".", 1)[0]+"_rare.vcf"

with open(gnomad_cohort_rare) as f2:
	ref_variants = f2.readlines()
	print(ref_variants)
	with open(disease_cohort) as f:
		for line in f:
			if line[0].isdigit():
'''				
