import  os, sys, subprocess, re

#SCRIPT FOR CHECKING HOW DISEASE IN CLINVAR
input_file = sys.argv[1]	#clinvar vcf

print(len(sys.argv), sys.argv)

if len(sys.argv) > 2:
	input_file2 = sys.argv[2]	#subset of genes to check for
	gene_list = []

	with open(input_file2, "r") as f:
		for line in f:
			line = line.strip()
			gene_list.append(line)
else:
	pass



#print(gene_list)

dis_list = []

with open(input_file, "r") as f:
	if len(sys.argv) > 2:
		for line in f:
			if "CLNDN=" in line and "GENEINFO=" in line:
				# print(line)
				gene = line.split("GENEINFO=")[1]
				gene = gene.split(":")[0]

				if gene in gene_list:

					disease = line.split("CLNDN=")[1]
					disease = disease.split(";")[0]
					disease = disease.split("|")
					for dis in disease:
						if not dis in dis_list:
							dis_list.append(dis)
	else:
		for line in f:
			if "CLNDN=" in line:
				disease = line.split("CLNDN=")[1]
				disease = disease.split(";")[0]
				disease = disease.split("|")
				for dis in disease:
					if not dis in dis_list:
						dis_list.append(dis)


if len(sys.argv) > 2:
	out_file = open(input_file+"_diseaselist_genesubset", "w")
else:
	out_file = open(input_file+"_diseaselist", "w")

out_file.write(str(len(dis_list))+"\n")
for x in dis_list:
	out_file.write(str(x)+"\n")
out_file.close() 

