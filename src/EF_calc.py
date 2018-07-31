import shlex, subprocess, os, sys, re

ID_list_file = sys.argv[1] #List of variants IDs that are present
CSV_file = sys.argv[2] #csv file e.g. HCM_missense_LMM_OMGL_chr.csv

ID_list = []
with open(ID_list_file) as f:
	for line in f:
		ID_list.append(line.rstrip())
print(ID_list)
cases_tot = 0
with open(CSV_file) as f:
	for line in f:
		line = line.split(",")
		if line[0] in ID_list:
			cases_tot += int(line[4])

print(cases_tot)