import  os, sys, subprocess, re

#SCRIPT FOR CHECKING HOW DISEASE IN CLINVAR
input_file = sys.argv[1]	#clinvar vcf

dis_list = []

with open(input_file, "r") as f:
	for line in f:
		if "CLNDN=" in line:
			line = line.split("CLNDN=")[1]
			line = line.split(";")[0]
			line = line.split("|")
			for dis in line:
				if not dis in dis_list:
					dis_list.append(dis)

out_file = open(input_file+"_diseaselist", "w")
out_file.write(str(len(dis_list)))
out_file.write(str(dis_list))
out_file.close() 