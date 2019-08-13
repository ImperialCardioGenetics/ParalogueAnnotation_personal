import  os, sys, subprocess, re,  shlex, subprocess, codecs, pysam, gzip, argparse

#SCRIPT FOR FILTERING OUT RARE VARIANTS DEFINED AS NOT BEING IN discovEHR

#to tabix index files - tabix -p vcf input.vcf.gz

input_file1 = sys.argv[1]   #path of tabix indexed file (.tab.gz file) #NOT tbi file
# input_file2 = sys.argv[2]  #path of query variants tableized input file 
tabix_file = pysam.TabixFile(input_file1)
# out_file_name = input_file1+"_w_REVEL"

# out_file = open(out_file_name, "w")


discovEHR_file = pysam.TabixFile("/work/nyl112/data/DiscovEHR/discovEHR_GRCh37.vcf.gz")

# def get_info_value(db_file,chrom,pos,ref_nt,alt_nt,ref_aa,alt_aa):
#     # position_found = False
#     REVEL= ""
#     for row in db_file.fetch(chrom, pos - 1, pos):
#         row_column = row.split('\t')
#         if str(pos) != row_column[1]:
#             continue
#         # position_found = True
#         row_ref_nt = row_column[2]
#         row_alt_nt = row_column[3]
#         row_ref_aa = row_column[4]
#         row_alt_aa = row_column[5]
#         #print row_ref_nt,row_alt_nt,row_ref_aa,row_alt_aa
#         #print ref_nt,alt_nt,ref_aa,alt_aa
#         if row_ref_nt == ref_nt and row_alt_nt == alt_nt and row_ref_aa == ref_aa and row_alt_aa ==alt_aa:
#             REVEL = row_column[6]
#         elif row_ref_nt == ref_nt and row_alt_nt == alt_nt:
#             REVEL = row_column[6]
#     # print(REVEL)
#     while not isinstance(REVEL, str):
#         REVEL = REVEL[0]
#     REVEL = [REVEL]
#     return REVEL

chrom = 1
pos = 69082


for pos in range(69080,69083):
	check = 0
	for row in tabix_file.fetch(chrom, pos - 1, pos):
		print(row)
		row_pos = int(row.split()[1])
		print(pos, row_pos)
		if pos == row_pos:
			print("TRUE")
			check = 1
		elif pos != row_pos:
			print("FALSE")

	if check == 1:
		print("CHECK PASSED")
	else:
		print("CHECK FAILED")

        # if line.startswith("CHROM"):
        # 	out_file.write(line.rstrip()+"\tREVEL_score\n")
        # if not line.startswith("CHROM"):
            # line = line.split()
            # if not "/" in line[8]:
            # 	out_file.write(str("\t".join(line+["NA"]))+"\n")
            # 	continue
            # chrom = line[0]
            # pos = int(line[1])
            # ntref = line[2]
            # ntalt = line[3]
            # aaref = line[8].split("/")[0]
            # aaalt = line[8].split("/")[1]
            # print(chrom,pos,ntref,ntalt,aaref,aaalt)
            # REVEL_score = get_info_value(tabix_file,chrom,pos,ntref,ntalt,aaref,aaalt)
            # print(REVEL_score)
            # new_line = "\t".join(line+REVEL_score)
            # out_file.write(str(new_line)+"\n")

# tabix_file.close()
# out_file.close()















'''
input_file = sys.argv[1]	#path of discovEHR vcf file to convert to vcf e.g. /media/nick/Data/PhD/Paralogues/Henrike's_data/discovEHR/GHS_Freeze_50.L3DP10.pVCF.frq.vcf
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

# try:
# 	gene = sys.argv[3]
# except IndexError:
# 	pass

with open(input_file) as f:
	ID_no = 1
	# header_out_file3_check = 0
	for og_line in f:
		# print(line)
		# if header_out_file3_check == 0:
		# 	out_file3.write("Gene,"+og_line)
		# 	header_out_file3_check = 1
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

			#if allele_freq <= 8.2e-06: #For LQTS according to https://www.nature.com/articles/gim201726/tables/1
			if allele_freq <= 1e-05: #For BrS

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
				# str(gene)+","+
				str(chrom)+","+
				str(pos)+","+
				str(ID)+","+
				str(csv_line[3])
				)

out_file.close()
out_file2.close()
out_file3.close()
'''















