def search_for_lines(filename, words_list):
    words_found = 0
    with open(filename) as db_file:
        for line_no, line in enumerate(db_file):
            if any(word in line for word in words_list):
                print(line_no, ':', line)
                words_found += 1
    return words_found

search_for_lines(
	"/media/nick/821ED5711ED55F2B/Users/Nick/Documents/PhD/Paralogues/ParalogueAnnotation_personal/data/clinvar/clinvar_20190114_GRCh38_onlyPathogenic_and_Likely_pathogenic.out_paraloc",
	
	)