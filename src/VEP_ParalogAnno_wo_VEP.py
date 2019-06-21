import os, sys, subprocess, re, shlex, subprocess, codecs


#Same as VEP_ParalogAnno.py but doesn't run VEP+Plugin as it assumes that has already been run. Saves time if .out_paraloc file has already by generated. See Paralogue_Annotation_masterscript_afterrun.py

def VEP_Plugin_afterrun(input_file, flavour=2, genome_build="GRCh38", VEPversion=90, offline=1, output_filename="default"):
	input_dir = input_file.rsplit("/", 1)[0]
	if flavour == 0:
		flavour = ""
		output_file = input_file.rsplit(".",1)[0]+".out_no_plugin"
		print(output_file)
	elif flavour == 1:
		flavour = "--plugin ParalogueAnnotation"
		output_file = input_file.rsplit(".",1)[0]+".out"
		print(output_file)
	elif flavour == 2:
		flavour = "--plugin ParalogueAnnotation,paraloc"
		output_file = input_file.rsplit(".",1)[0]+".out_paraloc"
		print(output_file)

	if output_filename != "default":
		output_file = output_filename	

	if not input_dir.endswith("/"):
		input_dir = input_dir+"/"
	print(input_dir)

	#Create custom IDs here
	out_file = open(
		output_file + "_tmp", "w", encoding="utf-8")
	with codecs.open(	#changed "open" to "codecs.open" to solve unicode error. Changed back maybe
		output_file, encoding="utf-8"
		) as infile:
		ID_no = 1
		for line in infile:
			if not line.startswith("#"):
				variant_info_line = line.split(None,3)
				variant_id = variant_info_line[2]
				if variant_id == ".":
					ID = "custom_" + str(ID_no)
					ID_no += 1
				else:
					ID = variant_id
				out_file.write(str(variant_info_line[0])+"\t"+str(variant_info_line[1])+"\t"+str(ID)+"\t"+str(variant_info_line[3]))
			else:
				out_file.write(line)
	out_file.close()

	#reads through outfile of vep+plugin and tidies up data by extracting only variants that have paralogous positions and spits that out as another outfile
	out_file = open(
		output_file + "_paralogs", "w", encoding="utf-8")

	output_file = output_file + "_tmp"

	if not flavour == "":
		with codecs.open(	#changed "open" to "codecs.open" to solve unicode error. Changed back maybe
			output_file, encoding="utf-8"
			) as infile:
			# ID_no = 1
			for line in infile:
				if not line.startswith("#"):
					line = line.split("CSQ=")
					# variant_info_line = line[0].split(None,3)
					CSQ = line[0]
					# print(CSQ)	#print query variant #commented out to avoid unicode error
					# variant_id = variant_info_line[2]
					# if variant_id == ".":
					# 	ID = "custom_" + str(ID_no)
					# 	ID_no += 1
					# else:
					# 	ID = variant_id
					line = line[1].split(",")
					# print(line, "\n")
					for x in line:
						y = x.split("|")
						gene = y[3]
						paralog = y[25]	#Remember to change indexing if also changing the number of flags being used by VEP
						if paralog == "" or paralog == "\n":
							paralog_check = 0
						else:
							paralog_check = 1
						if paralog_check == 1:
							# print(x, "\n")	#commented out to avoid unicode error
							out_file.write(
								str(CSQ)
								# str(variant_info_line[0])+"\t"+str(variant_info_line[1])+"\t"+str(ID)+"\t"+str(variant_info_line[3])
								+"PARALOGS->"
								+str(x)+"\n"
								)
				else:
					out_file.write(line)
	out_file.close()
	
	os.system("rm "+output_file)

	print("VEP_ParalogAnno Done!")
