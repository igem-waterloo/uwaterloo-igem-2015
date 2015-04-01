## Protocol followed to generate sgRNA sequences

1. Downloaded all *Caulimovrius* sequences from NCBI (includes CaMV).
   These sequences are in the file **caulimovirus_sequence.fasta**.
2. Ran through [Guidance v2](http://guidance.tau.ac.il/ver2/) to generate multiple sequence alignment.
3. Masked multiple sequence alignment based on base pairs that had at least 0.93% identity
   across all *Caulimovirus*.
4. Uploaded the masked sequences for CaMV to Benchling, which removed all gaps:
   [read-only link](https://benchling.com/s/If62yXXo/edit).
5. Ran Benchling CRISPR design, using parameters:
      * entire masked sequence set from selection
      * *A. thaliania* genome
      * Wild-Type Cas9 NGG PAM
      * 20 bp guide length
      * Genome region None
6. Kept sequences with efficiency score > 0.6 (calculated by [Benchling](https://benchling.com/) based on
   [Hsu et al.](http://crispr.mit.edu/about)) specificity score > 0.98 (calculated by
   Benchling based on [Doench et al.](http://www.nature.com/nbt/journal/v32/n12/abs/nbt.3026.html)).
7. Exported as primers and sanity-checked against positions in CaMV genome.
