---
title: "Paralog Annotation Notes"
# output: rmarkdown::github_document
output: 
  html_document:
    keep_md: true

editor_options: 
  chunk_output_type: inline
bibliography: bibliography.bib
always_allow_html: yes
---

<!--Load Packages and function-->


### Aims
* Apply Paralogue annotation on other datasets
    + "Genome Wide" - Clinvar dataset and ALL possible exome variants (get from `/data/Mirror/ExAC_release/release0.3.1/manuscript_data/all_possible_variants`)
    + Cardiomyopathy genes - [MYH7, MYBPC3, TNNT2, TPM1, MYL2, MYL3, TNNI3, ACTC1] -CHECK CODE TO SEE IF ALL MISSENSE CM VARIANTS REALLY DO NOT APPEAR IN ANY OF THESE GENES AS EG. MYH6 SHOULD HAVE SOME
    + Channelopathy genes - [KCNQ1, KCNH2, SCN5A, KCNE1, KCNE2, RYR2]
* Improve precision via increasing conservation of ref/alt alleles
    + pairwise QC - ignore any individual pairwise alignments where ref alleles are not conserved
    + pairwise QC and family QC - ignore entire alignment columns if the entire family ref alleles are not conserved; NB analogous to para z scores
    + pairwise QC and family QC and alt allele QC - ignore entire alignment columns if family ref allele and alt allele isn't conserved
* investigate para z scores
* investigate pfam meta domains

### Manuscript Plan
* New tool to show (more likely Erica will write up)
* Contrast to previous studies, is genome wide validated
* Describe implementation and how to use
    + vep plugin lib
* Provide additional descriptive statistics of input (clinvar) data, e.g. number of genes with paralogues, number of disease genes etc.
* Paralogue annotate P/LP with P/LP; paralogue annotate B/LB with P/LP
    + generate confusion matrices for above
* can PA also predict benign variants as well as pathogenic?
* does paraZ score add additional benefit
* additional test/validation dataset
    + disease
    + ExAC
* ICC genes - EFs 
* Distributability 
    + plugin
    + R shiny - vep web tool; integrated browser
    + integrated into gnomad
* Pfam domains

### Some interesting things to look at maybe
* "Ohnologs"
* Perform PA on paralogs from CPG (from @Modos2016) and non-CPG and see difference? 
* Integrating ortholog data to increase confidence calling 

### Introduction
With the advancements of sequencing technology, new potential variants are being discovered constantly. However to be able to identify said variants as pathogenic or benign requires supporting evidence, which does not always exists especially if the variant novel. 
Previously Ware *et al.* have developed **Paralogue Annotation** [@Ware2012; @Walsh2014], which utilizes information from paralogues (evolutionarily related genes from the same species) to help classify pathogenic variants. They verified its use in LQTS genes on variants acquired from patient cohorts.

Here **Paralogue Annotation** is tested further on a (Likely) Pathogenic/Benign varaint dataset from Clinvar.



### Material and Methods Notes

The Paralog Annotation algorithm was wriiten by Erica as a perl script plugin (called __ParalogueAnno_plugin_cleanup.pm__) for Ensembl's VEP version 90 (https://www.ensembl.org/info/docs/tools/vep/script/vep_options.html). 

The plugin has two arguments:

* the first parameter has 2 options:
    + ```variant``` (default) returns only the paralogous variants if any are present in the associated paralogs of the query gene found in the ensembl compara database
    + ```paraloc``` returns only paralog variant locations in the form of genomic coordinates of the corresponding codon in ALL paralogs;
* the second parameter has 2 options:
    + ```all``` for all variants;
    +  ```damaging``` (default) for only damaging variant.
The majority of the time ```paraloc``` mode is used.

The Ensembl team have touched up Erica's plugin and decrease runtime. The plugin is now called __ParalogueAnnotation.pm__

The initial output by VEP and the Plugin (VEP+Plugin) is not reader friendly for either the user nor if you want to parse informations. So a python wrapper, shwon below, for the VEP+Plugin was written to automatically parse the results, namely the paralogous variant information - __/data/Share/nick/Paralog_Anno/VEP_ParalogAnno.py__(Note the code below is not polished for release and is a WIP).


An intermediate python script (__File_prep_for_R.py__) was used to prep the results into R friendly data. Furthermore it could also be used to perform pairwise and family QC. Incidentally, the pairwise QC could be performed directly in R after the raw results are processed by __tableize_vcf.py__ and tabulated (see below).

__/data/Share/nick/Paralog_Anno/File_prep_for_R.py__ - formats results from __VEP_ParalogAnno.py__ into tabulated format ready for R processing


As paraloc mode only returns ref alleles. The alt alleles were extracted from the VEP information. This was done by using __tableize_vcf.py__.
__/data/Share/nick/Paralog_Anno/loftee/src/tableize_vcf.py__ was used to format the VEP output into table format for R processing. For example:

```bash
python /data/Share/nick/Paralog_Anno/loftee/src/tableize_vcf.py --vcf /data/Share/nick/Paralog_Anno/data_files/clinvar_20171029_onlyPathogenic.out_paraloc --out /data/Share/nick/Paralog_Anno/data_files/clinvar_20171029_onlyPathogenic.out_paraloc_tableized --do_not_minrep --include_id --vep_info SYMBOL,Amino_acids,Codons,Paralogue_Vars --split_by_transcript --canonical_only
```

If __--split_by_transcript__ is used then the code above is sufficient. Otherwise a python wrapper that includes additional formatting (__/data/Share/nick/Paralog_Anno/Tableize_wrapper.py__) that tableize couldn't do, i.e. separate variants that had multiple REF and ALT alleles was used to prepare the data for R. 


#### Datasets
Clinvar Likely Pathogenic/Pathogenic and Likely Benign/Benign variant vcf files were extracted and downloaded via the method developed by Zhang *et al.* [@Zhang2017]. 

#### Benchmarking performance of the plugin

__/data/Share/nick/Paralog_Anno/multi_vcf_extractor_benchmark.py__ is used to demonstrate speed at which VEP+Plugin takes to run


#### Scripts pipeline
<!-- old
vcf input file -> VEP_ParalogAnno.py -> File_prep_for_R.py -> Paralogous_var_align.R
vcf input file -> VEP_ParalogAnno.py -> paraloc file -> tableize_vcf.py (Tableize_wrapper.py) -> Paralogous_var_align.R
-->

```r
library(DiagrammeR)
pipeline = DiagrammeR::grViz("
digraph boxes_and_circles {
graph [overlap = true, fontsize = 10]

node [shape = plaintext, fillcolor = green, style=filled, fixedsize=false]
'VEP_ParalogAnno.py'; 'File_prep_for_R.py'; 'Tableize_wrapper.py'; 'R markdown'

node [shape = plaintext, fillcolor = orange, style=filled, fixedsize=false]
'vcf input file'; 'paralogs file'; 'paraloc file'; 'paralogs2 file'; 'paraloc_tableized file'

'vcf input file' -> 'VEP_ParalogAnno.py'; 'VEP_ParalogAnno.py' -> 'paralogs file'; 'VEP_ParalogAnno.py' -> 'paraloc file'; 'paralogs file' -> 'File_prep_for_R.py'; 'paraloc file' -> 'Tableize_wrapper.py'; 'File_prep_for_R.py' -> 'paralogs2 file'; 'Tableize_wrapper.py' -> 'paraloc_tableized file'; 'paralogs2 file' -> 'R markdown'; 'paraloc_tableized file' -> 'R markdown'

}")
pipeline
```

<!--html_preserve--><div id="htmlwidget-bb012f1c5951f59d0585" style="width:672px;height:480px;" class="grViz html-widget"></div>
<script type="application/json" data-for="htmlwidget-bb012f1c5951f59d0585">{"x":{"diagram":"\ndigraph boxes_and_circles {\ngraph [overlap = true, fontsize = 10]\n\nnode [shape = plaintext, fillcolor = green, style=filled, fixedsize=false]\n\"VEP_ParalogAnno.py\"; \"File_prep_for_R.py\"; \"Tableize_wrapper.py\"; \"R markdown\"\n\nnode [shape = plaintext, fillcolor = orange, style=filled, fixedsize=false]\n\"vcf input file\"; \"paralogs file\"; \"paraloc file\"; \"paralogs2 file\"; \"paraloc_tableized file\"\n\n\"vcf input file\" -> \"VEP_ParalogAnno.py\"; \"VEP_ParalogAnno.py\" -> \"paralogs file\"; \"VEP_ParalogAnno.py\" -> \"paraloc file\"; \"paralogs file\" -> \"File_prep_for_R.py\"; \"paraloc file\" -> \"Tableize_wrapper.py\"; \"File_prep_for_R.py\" -> \"paralogs2 file\"; \"Tableize_wrapper.py\" -> \"paraloc_tableized file\"; \"paralogs2 file\" -> \"R markdown\"; \"paraloc_tableized file\" -> \"R markdown\"\n\n}","config":{"engine":"dot","options":null}},"evals":[],"jsHooks":[]}</script><!--/html_preserve-->

#### Statistical terms
In context of is there a pathogenic paralogue alignment? A TP = pathogenic query variant with a paralogous pathogenic hit; FP = benign query variant with a paralogous pathogenic hit; FN = pathogenic query variant with no paralogous pathogenic hit; and TN= benign query variant with no paralogous pathogenic hit.

Likewise for a benign paralogous alignment, a TP = benign query variant with a paralogous benign hit; FP = pathogenic query variant with a paralogous benign hit; FN = benign query variant with no paralogous benign hit; and TN = pathogenic query variant with no paralogous benign hit.

#### Annotation of Clinvar
The Clinvar file __clinvar_20171029.vcf__ was downloaded from ftp://ftp.ncbi.nlm.nih.gov/pub/clinvar/vcf_GRCh38/. Note that since the initial look at what was available there's been updated Clinvar files. 

NOTE that I have noticed some descrepencies between the plugin annotations which call REFID = 1/0 and that of comparing the REF amino acid by VEP in the dataset to itself. This is due to the fact that the paralogous variant VEP is referring to is simply not in the dataset that I am annotating back to. As a result, it is best to make sure that the ref alleles are indeed the same when processing in R.

The annotataion of the entire clinvar set as of 2018:


<!--Previous old code below, remove once above code works-->
The annotataion of the entire clinvar set as of 20171029:


Taking only the 8 sarcomeric genes:

<!-- Old code below -->


Using only the 8 sarcomeric genes and joining to the whole clinvar dataset did not provide many annotations which could suggest either PA does not perform well on sarcomeric genes (paralogues to sarcomeric genes are not involed in disease) or that there is a lack of data. Therefore, it is not yet certain that PA does not work on sarcomeric genes and annotataion of additional sarcomeric data is required. See below.

Taking only the 5 channelopathy genes:


On the other hand, channelopathy genes did annotate well suggesting that their paralogues are involved in disease. 

Looking at alt alleles. Taking only pairwise alignments where the alt allele is conserved leaves only 1115 individual pairwise alignments. The number of actual unique variants this equates to is less - 825.

#### Annotation of all possible missense variations in the 8 sarcomeric genes and calculation of EF

For calculating the EFs, run the all possible missense variants through VEP+plugin and return paraloc locations. Then join those locations with pathogenic clinvar variants as before. This indicates which variants from all possible missense variants are likely to be pathogenic. Then we check to see if any of these variants are present in the cases and controls. Hopefully the controls will be less but there is more control data than cases bare in mind. Calculate the EFs using that. Remember though the EFs are based on how many times an allele is seen, not the number of different alleles by themselves.


<!--

-->

Total cases: 6140
number of affected cases: 39

Total controls: 60678
number of affected controls: 28

Odds ratio 	13.7648
95 % CI:	8.4648 to 22.3833
z statistic	10.570
Significance level	P < 0.0001
Attributable Risk Percent: 92.7% 
95 % CI: 79.6 to 100

#### Paralogue stats
According to ensembl, 92096 protein coding genes are defined to have paralogues. While 7958 protein genes do not have paralogues


#### Plan for overall annotations table
REF              | All | Alt matches |Alt no match
-----------------|-----|-------------|-------------
                 |     |             |
No QC            |     |             |
                 |     |             |
-----------------|-----|-------------|-------------
                 |     |             |
Paralog Conserved|     |             |
                 |     |             |
-----------------|-----|-------------|-------------
                 |     |             |
All Conserved    |     |             |
                 |     |             |
                 
#### Para-Z scores
For the para-z scores, will need to extract amino acid position from VEP output as well. Then look up the gene in question in para-z score folder, and using the position identify the para-z score. From my understanding, the para-z score is the same across aligned amino acids in the same gene family. Therefore, we could use a cut-off threshold to further improve our confidence in calling variants pathogenic etc. We could also then calculate ROC curves by altering the cut-off to see how that affects sensitivity/PPV.



#### Ohnologs
The "2R"" hypothesis states that some 500 million years ago, early vertebrates went through 2 rounds of whole genome duplication (WGD)[@Ohno1968]. Paralogues that arose from this WGD are known as ohnologs. @Singh2014 showed that monogenic disease genes to be enriched in ohnologs than other paralogs that arose from small scale duplications.
