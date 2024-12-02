# Project-Repository
This document contains all the necessary commands needed to replicate the analysis done for the project. All code containing gene family HADHA can be replaced with own gene family name.  

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
(Create new directory with gene family name, go to that directory, and make sure in that directory)

seqkit grep --pattern-file ~/lab03-$MYGIT/HADHA/HADHA.blastp.detail.filtered.out ~/lab03-$MYGIT/allprotein.fas | seqkit grep -v -p "carpio" > ~/lab04-$MYGIT/HADHA/HADHA.homologs.fas  
(Extract sequences in the "allprotein.fas" file that match the ones in the BLAST results file, and remove any sequences that contain "carpio", then save this to a new file)

less HADHA.homologs.fas  
(Use this to look at the newly created file to make sure the sequences were filtered correctly)

muscle -align ~/lab04-$MYGIT/HADHA/HADHA.homologs.fas -output ~/lab04-$MYGIT/HADHA/HADHA.homologs.al.fas  
(Align the sequences using MUSCLE and save the alignment in a new file)

alv -kli  ~/lab04-$MYGIT/HADHA/HADHA.homologs.al.fas | less -RS  
(View the alignment in a colored format, can see the differences very easily)

alv -kli --majority ~/lab04-$MYGIT/HADHA/HADHA.homologs.al.fas | less -RS  
(View the alignment, having the majority consensus sequence highlighted)

Rscript --vanilla ~/lab04-$MYGIT/plotMSA.R  ~/lab04-$MYGIT/HADHA/HADHA.homologs.al.fas  
(Make a visual plot of the alignment using R script)

alignbuddy -al  ~/lab04-$MYGIT/HADHA/HADHA.homologs.al.fas  
(This opens the alignment in the AlignBuddy tool for analysis, such as the wigth(length) of the alignment)

alignbuddy -trm all  ~/lab04-$MYGIT/HADHA/HADHA.homologs.al.fas | alignbuddy -al  
(Removes any gaps from the alignment, can view the length after this filtering)

alignbuddy -dinv 'ambig' ~/lab04-$MYGIT/HADHA/HADHA.homologs.al.fas | alignbuddy -al  
(This is used to calculate the length of the alignment after removing uncertain or ambiguous regions from the alignment)

t_coffee -other_pg seq_reformat -in ~/lab04-$MYGIT/HADHA/HADHA.homologs.al.fas -output sim  
(Using T-Coffee to calculate the average percent identity of the alignment)

alignbuddy -pi ~/lab04-$MYGIT/HADHA/HADHA.homologs.al.fas | awk ' (NR>2)  { for (i=2;i<=NF  ;i++){ sum+=$i;num++} }
END{ print(100*sum/num) }' 
(Using AlignBuddy to calculate the average percent identity of the alignment)

# Step 3: Creating Phylogenetic Tree 
Commands: 

mkdir ~/lab05-$MYGIT/HADHA
cd ~/lab05-$MYGIT/HADHA
pwd  
(Create new directory with gene family name, go to that directory, and make sure in that directory)

echo "(((((((((G.gallus_HADH_hydroxyacylcoenzyme_A_dehydrogenase_mitochondrial:0.2566511529, ...
...);" > ~/lab05-$MYGIT/HADHA.tree  
(Creates a newick-formatted species tree based on hardcorded branch lengths and topology, save it as a newly named file)

nw_display ~/lab05-$MYGIT/HADHA.tree  
(Displays the newick-formatted species tree of specified gene family, in this case HADHA)

nw_display -s ~/lab05-$MYGIT/HADHA.tree > ~/lab05-$MYGIT/HADHA.tree.svg 
(Generates graphic of the species tree and saves it into SVG format)

convert ~/lab05-$MYGIT/HADHA.tree.svg ~/lab05-$MYGIT/HADHA.tree.pdf 
(Converts the SVG species tree to a PDF)

Rscript --vanilla ~/lab05-$MYGIT/plotUnrooted.R ~/lab05-$MYGIT/HADHA.tree ~/lab05-$MYGIT/HADHA.unrooted.pdf 0.4 35  
(Using R script to plot the unrooted species tree based on the specified parameters and save it as a PDF file)

sed 's/ /_/g'  ~/lab04-$MYGIT/HADHA/HADHA.homologs.al.fas | seqkit grep -v -r -p "dupelabel" >  ~/lab05-$MYGIT/HADHA/HADHA.homologsf.al.fas  
(Remove sequences from the alignment that contains a duplicate tag or with the "dupelabel", replace the spaces with underscores in the alignment, save file as a copy)

iqtree -s ~/lab05-$MYGIT/HADHA/HADHA.homologsf.al.fas -bb 1000 -nt 2   
(Create a phylogenetic tree using IQ-Tree, containing 1000 bootstrap replicates and 2 CPU threads, this finds the maximum likehood tree estimate)

nw_display ~/lab05-$MYGIT/HADHA/HADHA.homologsf.al.fas.treefile  
(Displays the created phylogentic tree, newick formatted, view of the unrooted tree)

Rscript --vanilla ~/lab05-$MYGIT/plotUnrooted.R  ~/lab05-$MYGIT/HADHA/HADHA.homologsf.al.fas.treefile ~/lab05-$MYGIT/HADHA/HADHA.homologsf.al.fas.treefile.pdf 0.4 15 
(Generates a PDF of the unrooted tree using R script, allows the view of the tree with a graphical display, given specific scaling and margins based on the numerical values)

gotree reroot midpoint -i ~/lab05-$MYGIT/HADHA/HADHA.homologsf.al.fas.treefile -o ~/lab05-$MYGIT/HADHA/HADHA.homologsf.al.mid.treefile  
(Reroots the tree at the midpoint and saves this to a newly named file, midpoint rooting)

nw_order -c n ~/lab05-$MYGIT/HADHA/HADHA.homologsf.al.mid.treefile  | nw_display -  
(Reorders the nodes in the midpoint rooted tree, displays it so it can be seen, newick formatted)

nw_order -c n ~/lab05-$MYGIT/HADHA/HADHA.homologsf.al.mid.treefile | nw_display -w 1000 -b 'opacity:0' -s > ~/lab05-$MYGIT/HADHA/HADHA.homologsf.al.mid.treefile.svg -  
(Generate graphic of the tree, creating an SVG file of the rerooted tree with specified settings)

convert  ~/lab05-$MYGIT/HADHA/HADHA.homologsf.al.mid.treefile.svg  ~/lab05-$MYGIT/HADHA/HADHA.homologsf.al.mid.treefile.pdf  
(Converts the generated graphic of the rerooted tree from SVG to PDF)

nw_order -c n ~/lab05-$MYGIT/HADHA/HADHA.homologsf.al.mid.treefile | nw_topology - | nw_display -s -w 1000 > ~/lab05-$MYGIT/HADHA/HADHA.homologsf.al.midCl.treefile.svg -  
(Switches the view to a cladogram, provides a collapsed tree for easier visualization, save as SVG)

convert ~/lab05-$MYGIT/HADHA/HADHA.homologsf.al.midCl.treefile.svg ~/lab05-$MYGIT/HADHA/HADHA.homologsf.al.midCl.treefile.pdf  
(Converts this collapsed tree from SVG to PDF)

nw_reroot ~/lab05-$MYGIT/HADHA/HADHA.homologsf.al.fas.treefile H.sapiens_HBG1_hemoglobin_subunit_gamma1 H.sapiens_HBG2_hemoglobin_subunit_gamma2 H.sapiens_HBB_hemoglobin_subunit_beta H.sapiens_HBD_hemoglobin_subunit_delta > ~/lab05-$MYGIT/HADHA/HADHA.homologsf.outgroupbeta.treefile  
(Reroot the tree using outgroup rooting, using specific sequences as an outgroup based on common ancestry, save as a newly named file)

nw_order -c n ~/lab05-$MYGIT/HADHA/HADHA.homologsf.outgroupbeta.treefile | nw_topology - | nw_display -s -w 1000 > ~/lab05-$MYGIT/HADHA/HADHA.homologsf.outgroupbeta.treefile.svg -  
(Generate a graphic of the newly created outgroup rooted tree and save it as a SVG)

convert ~/lab05-$MYGIT/HADHA/HADHA.homologsf.outgroupbeta.treefile.svg ~/lab05-$MYGIT/HADHA/HADHA.homologsf.outgroupbeta.treefile.pdf  
(Convert the outgroup rooted tree from an SVG to a PDF)

# Step 4: Reconciliation of Species and Gene Trees
Commands: 

mkdir ~/lab06-$MYGIT/HADHA
cd ~/lab06-$MYGIT/HADHA
pwd  
(Create new directory with gene family name, go to that directory, and make sure in that directory)

cp ~/lab05-$MYGIT/HADHA/HADHA.homologsf.al.mid.treefile ~/lab06-$MYGIT/HADHA/HADHA.homologfs.al.mid.treefile  
(Copies the tree file of the midpoint rooted tree to this directory)

java -jar ~/tools/Notung-3.0_24-beta/Notung-3.0_24-beta.jar -s ~/lab05-$MYGIT/species.tre -g ~/lab06-$MYGIT/HADHA/HADHA.homologsf.al.mid.treefile --reconcile --speciestag prefix --savepng --events --outputdir ~/lab06-$MYGIT/HADHA/  
(Conducts the reconcilitation of species tree of the gene family HADHA with the gene tree based on the midpoint rooted tree using Notung, this generates PNG and event files)

nw_display ~/lab05-$MYGIT/species.tre  
(This displays the reconciled species tree)

grep NOTUNG-SPECIES-TREE ~/lab06-$MYGIT/HADHA/HADHA.homologs.al.mid.treefile.rec.ntg | sed -e "s/^\[&&NOTUNG-SPECIES-TREE//" -e "s/\]/;/" | nw_display -  
(This extracts the reconciled tree that was generated from Notung and displays it, showing the node names)

python2.7 ~/tools/recPhyloXML/python/NOTUNGtoRecPhyloXML.py -g ~/lab06-$MYGIT/HADHA/HADHA.homologsf.al.mid.treefile.rec.ntg --include.species  
(This converts the reconciled tree into RecPhyloXML format)

thirdkind -Iie -D 40 -f ~/lab06-$MYGIT/HADHA/HADHA.homologsf.al.mid.treefile.rec.ntg.xml -o ~/lab06-$MYGIT/HADHA/HADHA.homologsf.al.mid.treefile.rec.svg  
(Generates a graphic of the reconciled tree, into an SVG)

convert -density 150 ~/lab06-$MYGIT/HADHA/HADHA.homologsf.al.mid.treefile.rec.svg ~/lab06-$MYGIT/HADHA/HADHA.homologsf.al.mid.treefile.rec.pdf 
(Converts the SVG file of the reconciled tree to a PDF)

less HADHA.homologsf.al.mid.treefile.rec.events.txt  
(This allows the view of the duplications and losses in numerical format, can be used to determine the score(cost) of the reconciliation)

# Step 5: Protein Domain Analysis
Commands: 

mkdir ~/lab08-$MYGIT/HADHA && cd ~/lab08-$MYGIT/HADHA
pwd  
(Create new directory with gene family name, go to that directory, and make sure in that directory)

sed 's/*//' ~/lab04-$MYGIT/HADHA/HADHA.homologs.fas > ~/lab08-$MYGIT/HADHA/HADHA.homologs.fas  
(Make a copy of the raw unaligned sequences, removing the * or stop codon, save it to the new directory)

rpsblast -query ~/lab08-$MYGIT/HADHA/HADHA.homologs.fas -db ~/data/Pfam/Pfam -out ~/lab08-$MYGIT/HADHA/HADHA.rps-blast.out -outfmt "6 qseqid qlen qstart qend evalue stitle" -evalue .0000000001  
(Conducts RPS-BLAST on the homologous sequences with specific output format and e-value parameters, this uses the Pfam database for the RPS-BLAST)

rpsblast -query ~/lab08-$MYGIT/HADHA/HADHA.homologs.fas -db ~/data/Pfam/Pfam -out ~/lab08-$MYGIT/HADHA/HADHA.rps-blast.out -outfmt "6 qseqid qlen qstart qend evalue stitle" -evalue .1  
(Same RPS-BLAST conducted, but verifying the effect of the change in the e-value and why it should be lower)

cp ~/lab05-$MYGIT/HADHA/HADHA.homologsf.outgroupbeta.treefile ~/lab08-$MYGIT/HADHA  
(Copies the outgroup rooted tree that was previously created to this directory)

cp ~/lab05-$MYGIT/HADHA/HADHA.homologsf.al.mid.treefile ~/lab08-$MYGIT/HADHA  
(Copies the midpoint rooted tree that was previously created to this directory)

Rscript --vanilla ~/lab08-$MYGIT/plotTreeAndDomains.r ~/lab08-$MYGIT/HADHA/HADHA.homologsf.al.mid.treefile ~/lab08-$MYGIT/HADHA/HADHA.rps-blast.out ~/lab08-$MYGIT/HADHA/HADHA.tree.rps.pdf  
(This generates a phylogenetic tree that contains domain information on the side or domain annotations, this is done using R script and saves the file as a PDF)

mlr --inidx --ifs "\t" --opprint cat ~/lab08-$MYGIT/HADHA/HADHA.rps-blast.out | tail -n +2 | less -S  
(This uses Miller or mlr to view the output generated from the RPS-BLAST, it generates a nicely formatted and interactive view of the output)

cut -f 1 ~/lab08-$MYGIT/HADHA/HADHA.rps-blast.out | sort | uniq -c  
(Counting the unique query IDs in the RPS-BLAST output, counting the protiens that do not have any annotations)

cut -f 6 ~/lab08-$MYGIT/HADHA/HADHA.rps-blast.out | sort | uniq -c  
(Allows for the count of the most commonly found Pfam domain annotations)

awk '{a=$4-$3; print $1, "\t", a;}' ~/lab08-$MYGIT/HADHA/HADHA.rps-blast.out | sort -k2nr  
(Determines the alignment length for each query and sorts it by descending length, used to determine which protein has the longest annotated protein domain)

cut -f 1,5 -d $'\t' ~/lab08-$MYGIT/HADHA/HADHA.rps-blast.out  
(Pulls out the e-values for each of the protein domains and can be used to determine which domain annotation has the best or worst e-value)


