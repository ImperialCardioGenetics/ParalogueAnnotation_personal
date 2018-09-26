import  os, sys, subprocess, re

#SCRIPT FOR CHECKING WHAT GENES ARE PRESENT IN CLINVAR VCF FILE OF CHOICE (USE *out_paralocs FILE)

input_file = sys.argv[1]	#name of clinvar vcf file
dir1 = input_file.rsplit("/", 1)[0]
out_file = open(input_file+"_geneslist", "w")

gene_list = []
with open(input_file, "r") as f:
	for line in f:
		if line != "\n" and not line.startswith("#"):
			print(line)
			try:
				gene = line.split("GENEINFO=")[1].split(":")[0]
			except IndexError:
				continue
			print(gene)
			gene_list.append(gene)
	unq_list = set(gene_list)
	print(unq_list)
	for i in unq_list:
		out_file.write(i+"\n")

out_file.close()