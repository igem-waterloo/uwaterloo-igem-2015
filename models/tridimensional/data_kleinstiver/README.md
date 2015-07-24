The spreadsheet **Kleinstiver_mutants.xls** was created based on Extended Data Figure 2 of Kleinstiver et al. (2015),
"Engineered CRISPR-Cas9 nucleases with altered PAM specificities", Nature, [doi:10.1038/nature14592]( http://www.nature.com/nature/journal/vaop/ncurrent/full/nature14592.html). The figure gives the specific amino acid
substitutions in the PI domain of spyCas9 that lead to specificity for two novel PAMs, 'NGC' and 'NGA' (in contrast to
the standard 'NGG'). The CSVs **Kleinstiver_mutants_NGA.csv** and **Kleinstiver_mutants_NGC.csv** contain the same data
with less highlighting.

**kleinstiver_csv_to_dict** parses the CSVs into dictionaries which are annotated using the original PI domain sequence
(stored in **cas9_mutants.py** as PI_domain), the amino acid group (stored in **cas9_mutants.py** as aa_group) and the
secondary structure region where the mutant occurs (stored in **cas9_mutants.py** as PI_sec_structure). The output of
that script was then copied into **cas9_mutants.py** as kleinstiver_mutants.
