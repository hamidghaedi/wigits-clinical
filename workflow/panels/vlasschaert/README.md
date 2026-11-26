# Whitelist files of the Vlasschaert et al. CHIP filter

The five files beside this one are copied unchanged from `whitelist_filter_20230531.zip` in
https://github.com/briansha/Annovar_Whitelist_Filter_WDL at commit 8af2f67, the annotation and filtering
repository named in the author notes of Vlasschaert C, Mack T, Heimlich JB, et al., Blood
2023;141(18):2214-2223, doi:10.1182/blood.2022018825. They are distributed under the MIT licence in
`LICENSE`, copyright 2021 Charlie Condon and Brian Sharber.

| File | Content |
|---|---|
| Full_CHIP_gene_list_08262022.txt | 77 genes, the reported mutation rule and the RefSeq accession each rule is written against |
| CHIP_missense_vars_cv_08262022.txt | 1,435 gene and amino-acid change pairs |
| CHIP_nonsense_FS_vars_cv_08262022.txt | 47 genes in which nonsense and frameshift variants count |
| CHIP_splice_vars_cv_08262022.txt | 29 genes in which splice-site variants count |
| whitelist_filter_rscript_Aug26.R | the filter itself, kept as the reference for workflow/chip_curation.py |

`workflow/chip_curation.py` reads the four data files directly. Nothing in this directory is edited by
hand; to update, replace the files from a newer release of the repository.
