Paralog Annotation Notes
================

### Manuscript Plan

-   New tool to show (more likely Erica will write up)
-   Contrast to previous studies, is genome wide validated
-   Describe implementation and how to use
    -   vep plugin lib
-   Provide additional descriptive statistics of input (clinvar) data, e.g. number of genes with paralogues, number of disease genes etc.
-   Paralogue annotate P/LP with P/LP; paralogue annotate B/LB with P/LP
    -   generate confusion matrices for above
-   can PA also predict benign variants as well as pathogenic?
-   does paraZ score add additional benefit
-   additional test/validation dataset
    -   disease
    -   ExAC
-   ICC genes - EFs
-   Distributability
    -   plugin
    -   R shiny - vep web tool; integrated browser
    -   integrated into gnomad
-   Pfam domains

### Aims

-   Apply Paralogue annotation on other datasets
    -   "Genome Wide" - Clinvar dataset
    -   Cardiomyopathy genes - \[MYH7, MYBPC3, TNNT2, TPM1, MYL2, MYL3, TNNI3, ACTC1\]
    -   Channelopathy genes - \[KCNQ1, KCNH2, SCN5A, KCNE1, KCNE2, RYR2\]
-   Improve precision via increasing conservation of ref/alt alleles
    -   pairwise QC - ignore any individual pairwise alignments where ref alleles are not conserved
    -   pairwise QC and family QC - ignore entire alignment columns if the entire family ref alleles are not conserved; NB analogous to para z scores
    -   pairwise QC and family QC and alt allele QC - ignore entire alignment columns if family ref allele and alt allele isn't conserved
-   investigate para z scores
-   investigate pfam meta domains

### Material and Methods Notes

The Paralog Annotation algorithm was wriiten by Erica as a perl script plugin (called **ParalogueAnno\_plugin\_cleanup.pm**) for Ensembl's VEP version 90 (<https://www.ensembl.org/info/docs/tools/vep/script/vep_options.html>).

The plugin has two arguments:

-   the first parameter has 2 options:
    -   `variant` (default) returns only the paralogous variants if any are present in the associated paralogs of the query gene found in the ensembl compara database
    -   `paraloc` returns only paralog variant locations in the form of genomic coordinates of the corresponding codon in ALL paralogs;
-   the second parameter has 2 options:
    -   `all` for all variants;
    -   `damaging` (default) for only damaging variant. The majority of the time `paraloc` mode is used.

The Ensembl team have touched up Erica's plugin and decrease runtime. The plugin is now called **ParalogueAnnotation.pm**

The initial output by VEP and the Plugin (VEP+Plugin) is not reader friendly for either the user nor if you want to parse informations. So a python wrapper, shwon below, for the VEP+Plugin was written to automatically parse the results, namely the paralogous variant information - **/data/Share/nick/Paralog\_Anno/VEP\_ParalogAnno.py**(Note the code below is not polished for release and is a WIP).

An intermediate python script (**File\_prep\_for\_R.py**) was used to prep the results into R friendly data. Furthermore it could also be used to perform pairwise and family QC. Incidentally, the pairwise QC could be performed directly in R after the raw results are processed by **tableize\_vcf.py** and tabulated (see below).

**/data/Share/nick/Paralog\_Anno/File\_prep\_for\_R.py** - formats results from **VEP\_ParalogAnno.py** into tabulated format ready for R processing

As paraloc mode only returns ref alleles. The alt alleles were extracted from the VEP information. This was done by using **tableize\_vcf.py**. **/data/Share/nick/Paralog\_Anno/loftee/src/tableize\_vcf.py** was used to format the VEP output into table format for R processing. For example:

``` bash
python /data/Share/nick/Paralog_Anno/loftee/src/tableize_vcf.py --vcf /data/Share/nick/Paralog_Anno/data_files/clinvar_20171029_onlyPathogenic.out_paraloc --out /data/Share/nick/Paralog_Anno/data_files/clinvar_20171029_onlyPathogenic.out_paraloc_tableized --do_not_minrep --include_id --vep_info SYMBOL,Amino_acids,Codons,Paralogue_Vars
```

However, a python wrapper plus additional formatting (**/data/Share/nick/Paralog\_Anno/Tableize\_wrapper.py**) that tableize couldn't do, i.e. separate variants that had multiple REF and ALT alleles was used instead to prepare the data for R.

#### Benchmarking performance of the plugin

**/data/Share/nick/Paralog\_Anno/multi\_vcf\_extractor\_benchmark.py** is used to demonstrate speed at which VEP+Plugin takes to run

#### Scripts pipeline

vcf input file -&gt; VEP\_ParalogAnno.py -&gt; File\_prep\_for\_R.py -&gt; Paralogous\_var\_align.R vcf input file -&gt; VEP\_ParalogAnno.py -&gt; paraloc file -&gt; tableize\_vcf.py -&gt; Paralogous\_var\_align.R

#### Annotation of Clinvar

The Clinvar file **clinvar\_20171029.vcf** was downloaded from <ftp://ftp.ncbi.nlm.nih.gov/pub/clinvar/vcf_GRCh38/>. Note that since the initial look at what was available there's been updated Clinvar files.

NOTE that I have noticed some descrepencies between the plugin annotations which call REFID = 1/0 and that of comparing the REF amino acid by VEP in the dataset to itself. This is due to the fact that the paralogous variant VEP is referring to is simply not in the dataset that I am annotating back to. As a result, it is best to make sure that the ref alleles are indeed the same when processing in R.

The annotataion of the entire clinvar set as of 20171029:

Taking only the 8 sarcomeric genes:

Taking only the 5 channelopathy genes:

Looking at alt alleles. Taking only pairwise alignments where the alt allele is conserved leaves only 1115 individual pairwise alignments. The number of actual unique variants this equates to is less - 825.

#### Annotation of all possible missense variations in the 8 sarcomeric genes and calculation of EF

#### Paralogue stats

According to ensembl, 92096 protein coding genes are defined to have paralogues. While 7958 protein genes do not have paralogues

#### Plan for overall annotations table

| REF               | All   | Alt matches   | Alt no match  |
|-------------------|-------|---------------|---------------|
|                   |       |               |
| No QC             |       |               |               |
|                   |       |               |
| ----------------- | ----- | ------------- | ------------- |
|                   |       |               |
| Paralog Conserved |       |               |               |
|                   |       |               |
| ----------------- | ----- | ------------- | ------------- |
|                   |       |               |
| All Conserved     |       |               |               |
|                   |       |               |
