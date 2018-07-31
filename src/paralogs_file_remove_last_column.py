import shlex, subprocess, os, sys, re

input_file = sys.argv[1] #paralogs2 file

out_file_name = input_file+"woLastCol"
out_file = open(input_file+"woLastCol", "w")
with open(input_file) as f:
	for line in f:
		line = line.rstrip()
		out_file.write(str(line)+"\n")

out_file.close()
os.system("mv " + out_file_name + " " + input_file)