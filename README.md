# MyGenome
Analyses for ABT480/CS485G genome assembly 

## 1. Analyses of sequence quality 
The F1 and R1 sequence datasets were analyzed using FASTQC: 
```bash
ssh -Y cjea222@cjea222.cs.uky.edu
cd MyGenome/A16
fastqc &
```
Load F1 and R1 datasets into GUI interface. 
Take screen shots of output files:

[F1screenshot.png](/data/F1screenshot.png)
[F2screenshot.png](/data/F2screenshot.png)

## 2. Ran trimmomatic
```bash
java -jar sequences/trimmomatic-0.38.jar PE -threads 2 -phred33 -trimlog file.txt MyGenome/A16/UFVPY184_1.fq.gz MyGenome/A16/UFVPY184_2.fq.gz UFVPY184_1_paired.fastq UFVPY184_1_unpaired.fastq UFVPY184_2_paired.fastq UFVPY184_2_unpaired.fastq ILLUMINACLIP:MyGenome/A16/adaptors.fasta:2:30:10 SLIDINGWINDOW:20:20 MINLEN:100
```

## 3. Count number of forward reads remaining
```bash
grep -c \@A00261 UFVPY184_1_paired.fastq 
```

## 4. Count the total number of bases in both the forward and reverse files. 
```bash
grep -v -e \@A00216 -e \+  -e \F UFVPY184_1_paired.fastq | awk '{print length($0)}' | paste -sd+ | bc
grep -v -e \@A00216 -e \+  -e \F UFVPY184_2_paired.fastq | awk '{print length($0)}' | paste -sd+ | bc
```

The sum of the two command outputs above is the total number of bases in both the forward and reverse files. 

*screenshot of total number*

## 5. Analyses of trimmed sequence quality
```bash
ssh -Y cjea222@cjea222.cs.uky.edu
cd MyGenome/A16
fastqc &
```

Load F1 and R1 datasets into GUI interface. 
Take screen shots of output files:

[F1_trimmed_screenshot.png](/data/F1_trimmed_screenshot.png)
[F2_trimmed_screenshot.png](/data/F2_trimmed_screenshot.png)

## 6. Create a personal working directory in the MCC supercomputer. 
```bash
ssh cjea222@mcc.uky.edu
cd /project/farman_24cs485g/
mkdir cjea222
```

## 7. Transfer forward and reverse trimmed sequences to new directory. 
```bash
ssh cjea222@cjea222.cs.uky.edu
scp UFVPY184_1_paired.fastq cjea222@mcc.uky.edu:~/project/farman_24cs485g/cjea222
scp UFVPY184_2_paired.fastq cjea222@mcc.uky.edu:~/project/farman_24cs485g/cjea222
```

## 8. Use Velvet Advisor to identify a suitable k-mer value for assembling a bacterial genome. 
Launch online tool: https://dna.med.monash.edu/~torsten/velvet_advisor/
Fill in the relevant fields with sequence reads data. 
Estimate the genome size to be 5.5 Mb.
Use a target K-mer coverage of 20. 

[VelvetAdvisor.jpg](/data/VelvetAdvisor.jpg)

## 9. Run VelvetH and VelvetG with suggested K-mer value. 
```bash
ssh cjea222@mcc.uky.edu
cd /project/farman_s24cs485g/cjea222/
velveth UFVPY184_velvet 79 -shortPaired  -fastq.gz -separate UFVPY184_1_paired.fastq UFVPY184_2_paired.fastq
velvetg UFVPY184_velvet
```
Record genomze size, number of scaffolds, nd N50 value from Log file. 

*screenshot of results*

## 10. Run VelvetOptimiser using a range of K-mer values. 
```bash
ssh cjea222@mcc.uky.edu
cd /project/farman_s24cs485g/cjea222/
velvetoptimiser -s 121 -e 201 -x 10 -d UFVPY184_velvet_optimal -f '-shortPaired -fastq.gz  -separate UFVPY184_1_paired.fastq UFVPY184_2_paired.fastq' -t 1
```
Record genomze size, number of scaffolds, nd N50 value from Log file. 

*screenshot of results*

## 11. Rerun VelvetOptimiser using a narrower k-mer range and step size of 2 for best possible assembly of dataset. 
```bash
ssh cjea222@mcc.uky.edu
cd /project/farman_s24cs485g/cjea222/
velvetoptimiser -s 93 -e 119 -x 2 -d UFVPY184_velvet_optimal_1 -f '-shortPaired -fastq.gz  -separate UFVPY184_1_paired.fastq UFVPY184_2_paired.fastq' -t 1
```
Optimized asssembly with velvet hash value of 109. 
Record genomze size, number of scaffolds, nd N50 value from Log file. 

*screenshot of results*

## 12. Rename sequence headers to a standard format.
```bash
ssh cjea222@mcc.uky.edu
cd /project/farman_s24cs485g/
perl /SCRIPTS/SimpleFastaHeaders.pl cjea222/UFVPY184.fasta UFVPY184
```

## 13. Run Benchmarking Using Single-Copy Orthologs (BUSCO) on fully optimized assemblies. 
This is a metric used to assess genome assembly completeness. 
```bash
ssh cjea222@mcc.uky.edu
cd /project/farman_s24cs485g/
sbatch /project/farman_s24cs485g/SLURM_SCRIPTS/BuscoSingularity.sh MyGenome.fasta
```

*screenshot of BUSCO score*

## 14. Check Version of Blast on VM. 
```bash
ssh cjea222@cjea222.cs.uky.edu
blastn -version 
```

## 15. Copy MyGenome assembly (UFVPY184) from MCC to the blast directory on VM.
```bash
ssh cjea222@cjea222.cs.uky.edu
scp cjea222@mcc.uky.edu /project/farman_s24cs485g/cjea222/UFVPY184_nh.fasta blast/
```

## 16. BLAST UFVPY184 against sequences of transposons in the Magnaporthe genome.
```bash
ssh cjea222@cjea222.cs.uky.edu
cd blast
blastn -subject UFVPY184_nh.fasta -query MoRepeats.fasta -out MoRepeats.MyGenome.BLASTn(0-11) -evalue 1e-20 -outfmt (0-11)
```

## 17. Investigate aligned sequences. 
Find lengths of longest contigs in assembly: 
```bash
ssh cjea222@cjea222.cs.uky.edu
cd blast
scp AppliedBioinfo@10.163.183.71:~/Desktop/SequenceLengths.pl blast/
perl SequenceLengths.pl UFVPY184_nh.fasta | sort -k2n
```

Sort alignments based on chromosomal position: 
```bash
cd blast 
awk '$2 ~/contigXXX/' MoRepeats.UFVPY184.BLASTn6 | sort -k9n
```
or
```bash
cd blast 
grep <Repeat> MoRepeats.UFVPY184.BLASTn6 | awk \'$2 ~ /contigXXX/' | awk '$9 > (start position) && $9 < (end position)' | sort -k9n
```

Identify alignments involving full-length repeats: 
```bash
cd blast 
grep <Repeat> MoRepeats.UFVPY184BLASTn6 | awk '$4 >= 5638'
```

## 18. BLAST UFVPY184 against mitochondrial DNA and export list of aligned contigs. 
```bash
ssh cjea222@cjea222.cs.uky.edu
cd blast
blastn -query MoMitochondrion.fasta -subject UFVPY184_nh.fasta -evalue 1e-50 -max_target_seqs 20000 -outfmt '6 qseqid sseqid slen length qstart qend sstart send btop' -out MoMitochondrion.UFVPY184.BLAST
awk '$4/$3 > 0.9 {print $2 ",mitochondrion"}' MoMitochondrion.UFVPY184.BLAST > MyGenome_mitochondrion.csv
```

*list of aligned contigs*

## 19. BLAST UFVPY184 against repeat-masked version of the B71 reference genome. 
```bash
ssh cjea222@cjea222.cs.uky.edu
scp AppliedBioinfo@10.163.183.71:~/Desktop/B71v2sh_masked.fasta blast/ 
cd blast
blastn -query B71v2sh_masked.fasta -subject UFVPY184_nh.fasta -evalue 1e-50 -max_target_seqs 20000 -outfmt '6 qseqid sseqid qstart qend sstart send btop' -out B71v2sh.UFVPY184.BLAST
```

## 20. Identify SNP variants between UFVPY184 and B71v2sh genome. 
```bash
ssh cjea222@mcc.uky.edu
cd /project/farman_s24cs485g/cjea222/
scp cjea222@cjea222.cs.uky.edu:~/blast/B71v2sh.UFVPY184.BLAST
cd ..
cp SCRIPTs/CallVariants.sh cjea222/
cd cjea222
sbatch CallVariants.sh B71v2sh.UFVPY184.BLAST
```

## 21. Remove any contigs in UFVPY184 less than 200 base pairs in length.
```bash
ssh cjea222@mcc.uky.edu
cd /project/farman_s24cs485g/
cp SCRIPTs/CullShortContigs.pl cjea222/
cd cjea222
perl CullShortContigs.pl UFVPY184_nh.fasta
```

## 24. Use gene annotations from a reference genome to generate a set of training data for SNAP:
Download the B71ref2.fasta genome and B71Ref2_a0.3.gff3 annotation file from the Farman Mac Desktop to the current directory:
```bash
ssh cjea222@cjea222.cs.uky.edu
cd genes/snap
scp AppliedBioinfo@10.163.183.71:~/Desktop/B71ref2.fasta .
scp AppliedBioinfo@10.163.183.71:~/Desktop/B71Ref2_a0.3.gff3 .
```

Append the genome fasta sequence to the end of the gff3 file: 
```bash
echo '##FASTA' | cat B71Ref2_a0.3.gff3 - B71Ref2.fasta > B71Ref2.gff3
```

Convert the GFF file containing the reference data to a custom format (ZFF):
```bash
maker2zff B71Ref2.gff3
```

Extract the genome regions containing unique genes. Ask for up to 1000 base pairs of intergenic sequence on both sides of each gene:
```bash
fathom genome.ann genome.dna -categorize 1000
```

Extract the genome, transcript, and protein sequences from these genes. Keep 1000 base pairs of context and also flip genes that are on the reverse strand:
```bash
fathom uni.ann uni.dna -export 1000 -plus
```

Train the HMM:
```bash
forge export.ann export.dna
```

Condense all relative files into a single file for use with runs of SNAP:
```bash
hmm-assembler.pl Moryzae . > Moryzae.hmm
```

## 25. Running SNAP to search genomes for predicted genes. 
Run SNAP by giving the name of your parameter file and your FASTA file:
```bash
ssh cjea222@cjea222.cs.uky.edu
cd genes/snap
snap-hmm Moryzae.hmm UFVPY184_final.fasta > UFVPY184-snap.zff
```
Convert from default ZFF output format: 
```bash
snap-hmm Moryzae.hmm UFVPY184_final.fasta -gff > UFVPY184-snap.gff2
```

*include UFVPY184-snap.gff2 *

## 26. Running AUGUSTUS to search genomes for predicted genes. 
Magnaporthe grisea is very closely related to this species, so no need to retrain AUGUSTUS. Rather than specifying a parameter file explicitly, use the name of one of the included species:
```bash
ssh cjea222@cjea222.cs.uky.edu
cd genes/augustus
augustus --species=magnaporthe_grisea --gff3=on --singlestrand=true --progress=true ../snap/UFVPY184_final.fasta > UFVPY184-augustus.gff3
```

*include UFVPY184-augustus.gff3*

## 27. Combining evidence from SNAP and AUGUSTUS with MAKER. 
Create Maker configuration files:
```bash
ssh cjea222@cjea222.cs.uky.edu
cd genes/maker
maker -CTL
```

Open maker_opts.ctl with a text editor:
``bash
nano maker_opts.ctl
```

Find the lines for the following options and edit them as shown:
genome=/home/cjea222@cjea222/genes/snap/UFVPY184_final.fasta
model_org= must be set to blank
repeat_protein= must be set to blank
snaphmm=/home/cjea222@cjea222/genes/snap/Moryzae.hmm
augustus_species=magnaporthe_grisea
keep_preds=1
protein=/home/cjea222@cjea222/genes/maker/genbank/ncbi-protein-Magnaporthe_organism.fasta

Run Maker and log errors:
```bash
maker 2>&1 | tee maker.log
```

Merge all results from Maker into one GFF file: 
```bash
gff3_merge -d UFVPY184.maker.output/UFVPY184_master_datastore_index.log -o UFVPY184-annotations.gff
```

*include UFVPY184-annotations.gff* 
