
#SELECT gene_symbol, mut_id, vs_case_count, mut_exac_count, mut_coding, mut_protein, mut_effect, mut_exac_frequency, rep_patients
#FROM ICC_MUTATIONS.report, ICC_MUTATIONS.gene, ICC_MUTATIONS.mutation, ICC_MUTATIONS.variant_study
#WHERE rep_ref_id = 665 AND gene_id=mut_gene_id AND mut_id=rep_mut_id;

SELECT mut_chromosome, mut_chr_start_source, mut_id, mut_coding
FROM ICC_MUTATIONS.report, ICC_MUTATIONS.gene, ICC_MUTATIONS.mutation
WHERE rep_ref_id = 665 AND gene_id=mut_gene_id AND mut_id=rep_mut_id AND mut_effect = "missense";# AND mut_coding NOT LIKE "%del%";

SELECT gene_symbol, mut_id, mut_coding, mut_protein, mut_effect, rep_patients
FROM ICC_MUTATIONS.report, ICC_MUTATIONS.gene, ICC_MUTATIONS.mutation
WHERE rep_ref_id=665 AND gene_id=mut_gene_id AND mut_id=rep_mut_id AND mut_effect = "missense";
SELECT gene_symbol, mut_coding, mut_protein, mut_effect, mut_exac_frequency, rep_patients
from report, gene, mutation
where rep_ref_id=665 and gene_id=mut_gene_id and mut_id=rep_mut_id;


SELECT * FROM ICC_MUTATIONS.report, ICC_MUTATIONS.gene, ICC_MUTATIONS.mutation;


SELECT * FROM ICC_MUTATIONS.report, ICC_MUTATIONS.gene;
SELECT DISTINCT rep_disease FROM ICC_MUTATIONS.report;