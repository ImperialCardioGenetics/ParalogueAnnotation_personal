import os, sys, subprocess, re, shlex, subprocess, codecs


#WRAPPER FOR VEP(VERSION 90) and PLUGIN

def VEP_Plugin_run(input_file, flavour=2, genome_build="GRCh38", VEPversion=90, offline=1, output_filename="default"):
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


	if genome_build == "GRCh37":
		genome_build = "GRCh37 --port 3337"


	if VEPversion == 90:
		#VEP version 90
		mv_command = (
			"mv /data/Share/nick/Paralog_Anno/homo_sapiens_non_use/" + str(genome_build) + " /data/Share/nick/Paralog_Anno/homo_sapiens"
			)

		#OFFLINE MODE
		offline_command = (
			"perl -I /data/Share/nick/Paralog_Anno/paralogueAnnotator /data/Install/ensembl-vep/vep --force_overwrite --vcf --allele_number --canonical --offline --cache --dir_cache /data/Share/nick/Paralog_Anno/" +
			" --assembly " + str(genome_build) +
			" -i " + 
			input_file +
			" -o " + 
			output_file + 
			" " + flavour
			)

		#ONLINE MODE		
		online_command = (
			"perl -I /data/Share/nick/Paralog_Anno/paralogueAnnotator /data/Install/ensembl-vep/vep --force_overwrite --vcf --allele_number --canonical --database" +
			" --assembly " + str(genome_build) +
			" -i " + 
			input_file +
			" -o " + 
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
		print(online_command)
		os.system(online_command)
	elif offline == 1:
		print(offline_command)
		os.system(offline_command)
	
	with open(
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
		output_file + "_paralogs", "w", encoding="utf-8")
	if not flavour == "":
		with codecs.open(	#changed "open" to "codecs.open" to solve unicode error. Changed back maybe
			output_file, encoding="utf-8"
			) as infile:
		ID_no = 1
			for line in infile:
				if not line.startswith("#"):
					line = line.split("CSQ=")
					variant_info_line = line[0].split(None,3)
					# CSQ = line[0]
					# print(CSQ)	#print query variant #commented out to avoid unicode error
					# variant_info_line = CSQ.split()
					variant_id = variant_info_line[2]
					if variant_id == ".":
						ID = "custom_" + str(ID_no)
						ID_no += 1
					else:
						ID = variant_id
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
								#str(CSQ)
								str(variant_info_line[0])+"\t"+str(variant_info_line[1])+"\t"+str(ID)+"\t"+str(variant_info_line[3])
								+"PARALOGS->"
								+str(x)+"\n"
								)
				else:
					out_file.write(line)
	out_file.close()
	
	print("VEP_ParalogAnno Done!")
