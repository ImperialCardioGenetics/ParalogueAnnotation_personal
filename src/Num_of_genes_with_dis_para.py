import os, sys, subprocess, re

# script to find Number of disease protein coding genes that have at least 1 paralogue and at least 1 P/LP variant aligning to a non-gap position of another disease gene

input_file0 = "clinvar_20190114_GRCh38_onlyPathogenic_and_Likely_pathogenic.out_paraloc_tableized"
clinvar_P_LP_tableized_genes = []
with open(input_file0, "r") as f:
	for line in f:
		if not line.startswith("CHROM"):
			line = line.split()
			# print(line)
			if line[7] != "NA":
				gene = line[6]
				# print(gene)
				if not gene in clinvar_P_LP_tableized_genes:
					clinvar_P_LP_tableized_genes.append(gene)
		# inp = input('Enter to continue')
print(len(clinvar_P_LP_tableized_genes))
input_file = "clinvar_20190114_GRCh38_onlyPathogenic_and_Likely_pathogenic.out_paraloc_paralogs"

# num_of_genes_w_dis_para = 0
gene_check_list = []
with open(input_file, "r", encoding="utf-8") as f:
	for line in f:
		if line != "\n" and not line.startswith("#"):
			# chrom = line.split()[0]
			# position = line.split()[1]
			# ID = line.split()[2]
			# Ref = line.split()[3]
			# Alt = line.split()[4]
			
			line = line.rstrip()
			paralogs = line.split("PARALOGS->")[1]
			gene = paralogs.split("|")[3]
			# print(gene)
			# out_file.write(str(chrom)+" "+str(position)+"\t"+str(ID)+"\t"+str(gene)+"\t"+str(Ref)+"\t"+str(Alt)+"\t")

			paralogs = paralogs.split("|&")[1]
			paralogs = paralogs.split("&")
			paralogs.remove("")
			
			for paralog in paralogs:
				paralog_gene = paralog.split(":")[0]
				# print(paralog_gene)
				if paralog_gene in clinvar_P_LP_tableized_genes:
					gene_check_list.append(gene)

print(len(set(gene_check_list)))
# 1831 Number of disease protein coding genes that have at least 1 paralogue and at least 1 P/LP variant aligning to a non-gap position of another disease gene
