# Project-Repository
This document contains all the necessary commands needed to replicate the analysis done for the project.  

# Step 1: Finding Homologs 
Commands:

mkdir ~/lab03-$MYGIT/HADHA        
(Create the working directory for gene family name)

cd ~/lab03-$MYGIT/HADHA        
(Go to that directory, in this case HADHA)

pwd        
(Make sure in correct directory)

ncbi-acc-download -F fasta -m protein "NP_000173.2"        
(Download the query protein with the correct accession number, in this case for HADHA)

ls NP_000173.2.fa        
less NP_000173.2.fa        
(Check if the file exists and used to make sure the correct information is presented in the file)

blastp -db ../allprotein.fas -query NP_000173.2.fa -outfmt 0 -max_hsps 1 -out HADHA.blastp.typical.out        
(Conducts the BLAST search for a specific protein, in this case, HADHA)

blastp -db ../allprotein.fas -query NP_000173.2.fa -outfmt "6 sseqid pident length mismatch gapopen evalue bitscore pident stitle" -max_hsps 1 -out HADHA.blastp.detail.out        
(Conducts the BLAST search with a more detailed output)

less -S HADHA.blastp.detail.out        
(View the file that has the detailed output)

grep -c H.sapiens HADHA.blastp.detail.out        
(Counts the number of total human or H. sapiens hits based on the detailed output)  

awk '{if ($6< 1e-30)print $1 }' HADHA.blastp.detail.out > HADHA.blastp.detail.filtered.out        
(Filters out the homologs that have an e-value of less than 1e-30, so provides the homologs with high-confidence matches)

wc -l HADHA.blastp.detail.filtered.out        
(Counts the number of hits after filtering out the homologs)

grep -o -E "^[A-Z]\.[a-z]+" HADHA.blastp.detail.filtered.out | sort | uniq -c        
(Determine the number of paralogs found in each species based on the HADHA protein)

# Step 2: Making Alignments
Commands: 

mkdir ~/lab04-$MYGIT/HADHA        
cd ~/lab04-$MYGIT/HADHA        
pwd

seqkit grep --pattern-file ~/lab03-$MYGIT/HADHA/HADHA.blastp.detail.filtered.out ~/lab03-$MYGIT/allprotein.fas | seqkit grep -v -p "carpio" > ~/lab04-$MYGIT/HADHA/HADHA.homologs.fas

less HADHA.homologs.fas

muscle -align ~/lab04-$MYGIT/HADHA/HADHA.homologs.fas -output ~/lab04-$MYGIT/HADHA/HADHA.homologs.al.fas

alv -kli  ~/lab04-$MYGIT/HADHA/HADHA.homologs.al.fas | less -RS

alv -kli --majority ~/lab04-$MYGIT/HADHA/HADHA.homologs.al.fas | less -RS

Rscript --vanilla ~/lab04-$MYGIT/plotMSA.R  ~/lab04-$MYGIT/HADHA/HADHA.homologs.al.fas

alignbuddy -al  ~/lab04-$MYGIT/HADHA/HADHA.homologs.al.fas

alignbuddy -trm all  ~/lab04-$MYGIT/HADHA/HADHA.homologs.al.fas | alignbuddy -al

alignbuddy -dinv 'ambig' ~/lab04-$MYGIT/HADHA/HADHA.homologs.al.fas | alignbuddy -al

t_coffee -other_pg seq_reformat -in ~/lab04-$MYGIT/HADHA/HADHA.homologs.al.fas -output sim

alignbuddy -pi ~/lab04-$MYGIT/HADHA/HADHA.homologs.al.fas | awk ' (NR>2)  { for (i=2;i<=NF  ;i++){ sum+=$i;num++} }
END{ print(100*sum/num) }'

