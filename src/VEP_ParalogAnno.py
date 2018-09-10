import  os, sys, subprocess, re, shlex, subprocess, codecs


#WRAPPER FOR VEP(VERSION 90) and PLUGIN

# def VEP_Plugin_run(flavour, input_dir, input_file, genome_build, VEPversion):
def VEP_Plugin_run(input_file, flavour=2, genome_build="GRCh38", VEPversion=90, offline=1, output_filename="default"):
	# flavour = int(sys.argv[1])	#0(no plugin), 1(variants), 2(paraloc)
	# input_dir = sys.argv[2]	#input directory
	input_dir = input_file.rsplit("/", 1)[0]
	# input_file = sys.argv[3]	#name of input file
	# # output_file = sys.argv[3]	#name of output file
	# genome_build = str(sys.argv[4]) #37, 38
	# VEPversion = int(sys.argv[5]) #83, 87, 90

	#e.g. python3 VEP_ParalogAnno.py 0 /data/Share/nick/Paralog_Anno/data_files 10_pathogenic.vcf 38 90

	# if not len(os.listdir("/data/Share/nick/Paralog_Anno/homo_sapiens/")) == 0:
	# 	sys.exit("Terminating script:/data/Share/nick/Paralog_Anno/homo_sapiens/ not empty. Please empty and try again")

	# print(type(flavour))
	if flavour == 0:
		# print("check")
		flavour = ""
		output_file = input_file.rsplit(".",1)[0]+".out_no_plugin"
		print(output_file)
	# elif flavour == 1:
	# 	flavour = "--plugin ParalogueAnno_plugin_cleanup"
	# 	output_file = input_file.rsplit(".",1)[0]+".out_old"
	# 	print(output_file)
	# elif flavour == 2:
	# 	flavour = "--plugin ParalogueAnno_plugin_cleanup,paraloc"
	# 	output_file = input_file.rsplit(".",1)[0]+".out_old_paraloc"
	# 	print(output_file)
	# elif flavour == 3:
	# 	flavour = "--plugin ParalogueAnnotation"
	# 	output_file = input_file.rsplit(".",1)[0]+".out_new"
	# 	print(output_file)
	# elif flavour == 4:
	# 	flavour = "--plugin ParalogueAnnotation,paraloc"
	# 	output_file = input_file.rsplit(".",1)[0]+".out_new_paraloc"
	# 	print(output_file)
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


	'''if VEPversion == 83:
		#VEP version 83
		# rm_command = ("rm -rf /data/Share/nick/Paralog_Anno/homo_sapiens/83_GRCh" + genome_build)
		mv_command = (
			"mv /data/Share/nick/Paralog_Anno/homo_sapiens_non_use/83_GRCh" + genome_build + " /data/Share/nick/Paralog_Anno/homo_sapiens"
			)
		command = (
			"perl -I /data/Install/ensembl-tools-release-83/scripts/variant_effect_predictor/.vep/Plugins /data/Install/ensembl-tools-release-83/scripts/variant_effect_predictor/variant_effect_predictor.pl --force_overwrite --vcf --offline --cache --dir_cache /data/Share/nick/Paralog_Anno/" +
			" -i " + input_dir + input_file +
			" -o " + input_dir + output_file + 
			" " + flavour
			)
		mv_back_command = (
			"mv /data/Share/nick/Paralog_Anno/homo_sapiens/83_GRCh" + genome_build + " /data/Share/nick/Paralog_Anno/homo_sapiens_non_use"
			)
	elif VEPversion == 87:
		#VEP version 87
		# rm_command = ("rm -rf /data/Share/nick/Paralog_Anno/homo_sapiens/87_GRCh" + genome_build)
		mv_command = (
			"mv /data/Share/nick/Paralog_Anno/homo_sapiens_non_use/87_GRCh" + genome_build + " /data/Share/nick/Paralog_Anno/homo_sapiens"
			)
		command = (
			"perl -I /data/Install/ensembl-tools-release-83/scripts/variant_effect_predictor/.vep/Plugins /data/Install/ensembl-tools-release-87/ensembl-tools/scripts/variant_effect_predictor/variant_effect_predictor.pl --force_overwrite --vcf --offline --cache --dir_cache /data/Share/nick/Paralog_Anno/" +
			" -i " + input_dir + input_file +
			" -o " + input_dir + output_file + 
			" " + flavour
			)
		mv_back_command = (
			"mv /data/Share/nick/Paralog_Anno/homo_sapiens/87_GRCh" + genome_build + " /data/Share/nick/Paralog_Anno/homo_sapiens_non_use"
			)'''
	if VEPversion == 90:
		#VEP version 90
		# rm_command = ("rm -rf /data/Share/nick/Paralog_Anno/homo_sapiens/90_GRCh" + genome_build)
		mv_command = (
			"mv /data/Share/nick/Paralog_Anno/homo_sapiens_non_use/" + str(genome_build) + " /data/Share/nick/Paralog_Anno/homo_sapiens"
			)

		#OFFLINE MODE
		offline_command = (
			# "perl -I /data/Install/ensembl-tools-release-83/scripts/variant_effect_predictor/.vep/Plugins /data/Install/ensembl-vep/vep --force_overwrite --vcf --offline --cache --dir_cache /data/Share/nick/Paralog_Anno/" +
			"perl -I /data/Share/nick/Paralog_Anno/paralogueAnnotator /data/Install/ensembl-vep/vep --force_overwrite --vcf --allele_number --canonical --offline --cache --dir_cache /data/Share/nick/Paralog_Anno/" +
			" --assembly " + str(genome_build) +
			" -i " + 
			# input_dir + 
			input_file +
			" -o " + 
			# input_dir + 
			output_file + 
			" " + flavour
			)

		#ONLINE MODE		
		online_command = (
			# "perl -I /data/Install/ensembl-tools-release-83/scripts/variant_effect_predictor/.vep/Plugins /data/Install/ensembl-vep/vep --force_overwrite --vcf --offline --cache --dir_cache /data/Share/nick/Paralog_Anno/" +
			"perl -I /data/Share/nick/Paralog_Anno/paralogueAnnotator /data/Install/ensembl-vep/vep --force_overwrite --vcf --allele_number --canonical --database" +
			" --assembly " + str(genome_build) +
			" -i " + 
			# input_dir + 
			input_file +
			" -o " + 
			# input_dir + 
			output_file + 
			" " + flavour
			)
		
		mv_back_command = (
			"mv /data/Share/nick/Paralog_Anno/homo_sapiens/" + str(genome_build) + " /data/Share/nick/Paralog_Anno/homo_sapiens_non_use"
			)

	elif VEPversion == 93:
		#OFFLINE MODE
		offline_command = (
			"perl -I /data/Share/nick/Paralog_Anno/paralogueAnnotator /data/Share/nick/Paralog_Anno/ensembl-vep/vep --force_overwrite --vcf --allele_number --canonical --offline --cache --dir_cache /data/Share/nick/Paralog_Anno/ --assembly " + 
			str(genome_build) + " -i " + input_file + " -o " + output_file + " " + flavour
			)
		#ONLINE MODE
		online_command = (
			"perl -I /data/Share/nick/Paralog_Anno/paralogueAnnotator /data/Share/nick/Paralog_Anno/ensembl-vep/vep --force_overwrite --vcf --allele_number --canonical --database --assembly " +
			str(genome_build) + " -i " + input_file + " -o " + output_file + " " + flavour			
			)
	
	if offline == 0:
		# print(mv_command)
		# os.system(mv_command)
		print(online_command)
		os.system(online_command)
		# print(mv_back_command)
		# os.system(mv_back_command)
	elif offline == 1:
		print(offline_command)
		os.system(offline_command)
	
	with open(
		# input_dir + 
		output_file + 
		"_summary.html"
		) as infile:
		for line in infile:
			if "<td>Run time</td>" in line:
				if VEPversion == 83 or VEPversion == 87:
					#IF VEP87 AND 83
					m = re.search(r"<td>Run time<\/td> <td>.*second.<\/td>", line)
					# print(m.group())
					run_time = m.group().split("</td> <td>")[1].rstrip("</td>")
					print(run_time, "\n")
				elif VEPversion == 90:
					# #IF VEP90
					m = re.search(r"<td>Run time<\/td><td>.*second.<\/td>", line)
					# print(m.group())
					run_time = m.group().split("</td><td>")[1].rstrip("</td>")
					print(run_time, "\n")

	###TESTED THE BELOW WITH VEP90#### MIGHT NOT WORK WITH OTHER VERSIONS
	#reads through outfile of vep+plugin and tidies up data by extracting only variants that have paralogs and spits that out as another outfile
	out_file = open(
		# input_dir + 
		output_file + "_paralogs", "w", encoding="utf-8")
	if not flavour == "":
		with codecs.open(	#changed "open" to "codecs.open" to solve unicode error. Changed back maybe
			# input_dir + 
			output_file, encoding="utf-8"
			) as infile:
			for line in infile:
				if not line.startswith("#"):
					line = line.split("CSQ=")
					CSQ = line[0]
					# print(CSQ)	#print query variant #commented out to avoid unicode error
					line = line[1].split(",")
					# print(line, "\n")
					for x in line:
						y = x.split("|")
						# print(y, "\n")
						gene = y[3]
						paralog = y[25]	#Remember to change indexing if also changing the number of flags being used by VEP
						if paralog == "" or paralog == "\n":
							paralog_check = 0
						else:
							paralog_check = 1
						if paralog_check == 1:
							# print(x, "\n")	#commented out to avoid unicode error
							out_file.write(
								str(CSQ)+"PARALOGS->"
								+str(x)+"\n"
								)
				else:
					out_file.write(line)
	out_file.close()
	
	# ###FOR CLINVAR VARIANTS AT LEAST###SCRAP ALL OF THIS, USE R TIDYVERSE JOIN METHOD INSTEAD
	# out_file = open(input_dir + output_file + "_aligned_vars", "w")
	# with open(input_dir + output_file + "_paralogs") as infile1:
	# 	infile2 = open(input_dir + input_file, "r") 
	# 	text_list = infile2.read()
	# 	infile2.close()
	# 	####TEST_BELOW####
	# 	# chromosome = 1
	# 	# pos1 = 1050574
	# 	# pos2 = 1050575
	# 	# pos3 = 1050576
	# 	# s1 = re.findall("^"+str(chromosome)+"\t"+str(pos1)+"\s.*$", text_list, re.MULTILINE)
	# 	# s2 = re.findall("^"+str(chromosome)+"\t"+str(pos2)+"\s.*$", text_list, re.MULTILINE)
	# 	# s3 = re.findall("^"+str(chromosome)+"\t"+str(pos3)+"\s.*$", text_list, re.MULTILINE)		
	# 	# if s1:
	# 	# 	print("YES")
	# 	# 	print(s1)
	# 	# if s2:
	# 	# 	print(s2)
	# 	# if s3:
	# 	# 	print(s3)
	# 	####END_OF_TEST####
	# 	for line in infile1:
	# 		if not line.strip() == "":
	# 			query_variant = line.split("PARALOGS->")[0]
	# 			out_file.write(query_variant)
	# 			paralog_line = line.split("PARALOGS->")[1].rstrip()
	# 			paralog_line = paralog_line.split("|&")[1]
	# 			paralog_line = paralog_line.split("&")
	# 			paralog_line = list(filter(None, paralog_line))
	# 			print(paralog_line)
	# 			out_file.write(str(paralog_line)+";")
	# 			for paralog in paralog_line:
	# 				if not paralog == "":
	# 					print(paralog)
	# 					location = paralog.split(":")[1]
	# 					chromosome = location.split("_")[0].lstrip("chr")
	# 					position = location.split("_")[1]
	# 					position1 = int(position.split("-")[0])
	# 					position2 = int(position.split("-")[1])
	# 					positions = list(range(position1,position2+1))
	# 					print(chromosome, type(positions), positions)
	# 					# p = re.compile(r"^"+str(chromosome)+r"\t"+positions+r"\t")
	# 					for pos in positions:
	# 						s = re.findall("^"+str(chromosome)+"\t"+str(pos)+"\s.*$",text_list, re.MULTILINE)
	# 						if s:
	# 							print(s)
	# 							out_file.write(str(s)+";")
	# 			out_file.write("\n")
	# out_file.close()
	# #################################################
	
	print("VEP_ParalogAnno Done!")

# VEP_Plugin_run(2, "/data/Share/nick/Paralog_Anno/data_files", "new_extracted_clinvar_cm_intepretated_raw_variants_2.vcf", "38", 87)


# flavour = int(sys.argv[1])	#0(no plugin), 1(variants), 2(paraloc),
# input_dir = sys.argv[2]	#input directory
# input_file = sys.argv[3]	#name of input file
# genome_build = str(sys.argv[4]) #37, 38
# VEPversion = int(sys.argv[5]) #83, 87, 90
# VEP_Plugin_run(flavour, input_dir, input_file, genome_build, VEPversion)
