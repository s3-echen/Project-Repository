# Project-Repository
This document contains all the necessary commands needed to replicate the analysis done for the project.  

# Step 1: Finding Homologs 
Commands:
mkdir ~/lab03-$MYGIT/HADHA
//Create the working directory for gene family name
cd ~/lab03-$MYGIT/HADHA
//Go to that directory, in this case HADHA
pwd
//Make sure in correct directory
ncbi-acc-download -F fasta -m protein "NP_000173.2"
ls NP_000173.2.fa
less NP_000173.2.fa
blastp -db ../allprotein.fas -query NP_000173.2.fa -outfmt 0 -max_hsps 1 -out HADHA.blastp.typical.out
blastp -db ../allprotein.fas -query NP_000173.2.fa -outfmt "6 sseqid pident length mismatch gapopen evalue bitscore pident stitle" -max_hsps 1 -out HADHA.blastp.detail.out
less -S HADHA.blastp.detail.out
grep -c H.sapiens HADHA.blastp.detail.out
awk '{if ($6< 1e-30)print $1 }' HADHA.blastp.detail.out > HADHA.blastp.detail.filtered.out
wc -l HADHA.blastp.detail.filtered.out
grep -o -E "^[A-Z]\.[a-z]+" HADHA.blastp.detail.filtered.out | sort | uniq -c

