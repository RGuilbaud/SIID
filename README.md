# Scripts for : Higher levels of SARS-CoV-2 genetic variation among immunocompromised patients: a retrospective case-control study

Mapping and pileup generation 
  -	Mapping, filtering, formatting and creation of pileup file then tsv: script_bamcov_gz.sh + extract_info_from_pileup.pl 
  -	Extract number of fasta lines before and after filtering: read_count.sh (! not actual read count but number of lines: divide these values by 4)

Mutation detection and analysis 
  -	Coverage analysis:  coverage per sample + Ct graph, nb reads and median coverage: 0_sample_coverage _bichat_accurate_seq1-7.R  + 0_sample_coverage_pitie.R  
  -	Split tsv by gene: 1_tsv_per_gene_bichat_accurate.R
  -	Search for mutations (by gene) and concatenate results: 2_finding_mutations_RG_V2_bichat_accurate.R + Fonction_code_genetique.R and   3_Mutations_recap_table_bichat_accurate.R
  -	Create a multifasta of consensus sequences: 4_Make_fasta_from_mutation_file.R + 4_Make_fasta_from_mutation_file_pitie.R
  -	Obtain variants from samples: Nextclade 
  -	Creation of sample summary tables + metadata: 4_Infos_tab_seq1234567.R
  -	Creation of a file containing sample info for each mutation (genot matrix equivalent): 4_2_Mutation_table_full_bichat_accurate.R 
  -	Change the major/minor threshold to 20% (it is 50% for consensus): 4_3_change_limit_minor_major.R
  -	Creation of graphs describing mutations (including heatmap and PCA) 
    •	5_graphes_mutations_Bichat_accurate.R 
    •	Missing data filter + PCA and heatmap: 5_Missing_data_PCA_heatmap.R
  - Deletion analysis: 5Bis_del_analysis.R
  -	Position of mutations on the genome: 6_Position_graph_mutations_full_genome.R (complete genome) 
  -	Patient vs. staff comparison: 7_Comp_patients_personnel.R 
  -	Follow-up analysis : 8_Patients_suivis.R
  -	Analysis of potential co-infections: 9_Potential_coinfections.R
