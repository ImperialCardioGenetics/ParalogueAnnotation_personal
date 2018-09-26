import  os, sys, subprocess, re
	
#QUICK SCRIPT FOR FINDING UNIQUE ITEMS (GENES) IN A LIST (GENELIST)

input_file = sys.argv[1]	#name of genelist file, expects to be a one column list of genes with every new gene on new line
dir1 = input_file.rsplit("/", 1)[0]
out_file = open(input_file+"_Unique", "w")

genelist = []
with open(input_file, "r") as f:
	for gene in f:
		gene = gene.rstrip()
		genelist.append(gene)


unq_list = set(genelist)
print(unq_list)
for gene in unq_list:
	out_file.write(gene+"\n")

out_file.close()