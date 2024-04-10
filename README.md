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

## 8. Copy Velvetoptimiser script to new directory. 
```bash
ssh cjea222@cjea222.cs.uky.edu
- use nano to edit

```
## 9. Use Velvet Advisor. 

## 10. Run VelvetH and VelvetG with suggested K-mer value. 
```bash
ssh cjea222@mcc.uky.edu
cd /project/farman_s24cs485g/cjea222/
velveth UFVPY184_velvet 79 -shortPaired  -fastq.gz -separate UFVPY184_1_paired.fastq UFVPY184_2_paired.fastq
velvetg UFVPY184_velvet
```
Record genomze size, number of scaffolds, nd N50 value from Log file. 

## 11. Run VelvetOptimiser using a range of K-mer values. 
```bash
ssh cjea222@mcc.uky.edu
cd /project/farman_s24cs485g/cjea222/
velvetoptimiser -s 121 -e 201 -x 10 -d UFVPY184_velvet_optimal -f '-shortPaired -fastq.gz  -separate UFVPY184_1_paired.fastq UFVPY184_2_paired.fastq' -t 1
```
Record genomze size, number of scaffolds, nd N50 value from Log file. 

## 12. Rerun VelvetOptimiser using a narrower k-mer range and step size of 2 for best possible assembly of dataset. 
```bash
ssh cjea222@mcc.uky.edu
cd /project/farman_s24cs485g/cjea222/
velvetoptimiser -s 93 -e 119 -x 2 -d UFVPY184_velvet_optimal_1 -f '-shortPaired -fastq.gz  -separate UFVPY184_1_paired.fastq UFVPY184_2_paired.fastq' -t 1
```
Optimized asssembly with velvet hash value of 109. 
Record genomze size, number of scaffolds, nd N50 value from Log file. 

## 13. Rename sequence headers to a standard format.
```bash
ssh cjea222@mcc.uky.edu
cd /project/farman_s24cs485g/
perl /SCRIPTS/SimpleFastaHeaders.pl cjea222/UFVPY184.fasta UFVPY184
```

## 14. Run Benchmarking Using Single-Copy Orthologs (BUSCO) on fully optimized assemblies. 
```bash
ssh cjea222@mcc.uky.edu
cd /project/farman_s24cs485g/
sbatch /project/farman_s24cs485g/SLURM_SCRIPTS/BuscoSingularity.sh MyGenome.fasta
```

## 15. Check Version of Blast on VM. 
```bash
ssh cjea222@cjea222.cs.uky.edu
blastn -version 
```

## 16. Copy MyGenome assembly (UFVPY184) from MCC to the blast directory on VM.
```bash
ssh cjea222@cjea222.cs.uky.edu
scp cjea222@mcc.uky.edu /project/farman_s24cs485g/cjea222/UFVPY184_nh.fasta blast/
```

## 17. BLAST UFVPY184 against sequences of transposons in the Magnaporthe genome.
```bash
ssh cjea222@cjea222.cs.uky.edu
cd blast
blastn -subject UFVPY184_nh.fasta -query MoRepeats.fasta -out MoRepeats.MyGenome.BLASTn(0-11) -evalue 1e-20 -outfmt (0-11)
```

## 18. Investigate aligned sequences. 
Find lengths of longest contigs in assembly: 
```bash
ssh cjea222@cjea222.cs.uky.edu
cd blast
scp AppliedBioinfo@10.163.183.71:~/Desktop/SequenceLengths.pl blast/
perl SequenceLengths.pl UFVPY184_nh.fasta | sort -k2n
```

Sort alignments based on chromosomal position: 
```bash
ssh cjea222@cjea222.cs.uky.edu
cd blast 
awk '$2 ~/contigXXX/' MoRepeats.UFVPY184.BLASTn6 | sort -k9n
```
or
```bash
ssh cjea222@cjea222.cs.uky.edu
cd blast 
grep <Repeat> MoRepeats.UFVPY184.BLASTn6 | awk \'$2 ~ /contigXXX/' | awk '$9 > (start position) && $9 < (end position)' | sort -k9n
```

Identify alignments involving full-length repeats: 
```bash
ssh cjea222@cjea222.cs.uky.edu
cd blast 
grep <Repeat> MoRepeats.UFVPY184BLASTn6 | awk '$4 >= 5638'
```

## 19. BLAST UFVPY184 against mitochondrial DNA and export list of aligned contigs. 
```bash
ssh cjea222@cjea222.cs.uky.edu
cd blast
blastn -query MoMitochondrion.fasta -subject UFVPY184_nh.fasta -evalue 1e-50 -max_target_seqs 20000 -outfmt '6 qseqid sseqid slen length qstart qend sstart send btop' -out MoMitochondrion.UFVPY184.BLAST
awk '$4/$3 > 0.9 {print $2 ",mitochondrion"}' MoMitochondrion.UFVPY184.BLAST > MyGenome_mitochondrion.csv
```

## 20. BLAST UFVPY184 against repeat-masked version of the B71 reference genome. 
```bash
ssh cjea222@cjea222.cs.uky.edu
scp AppliedBioinfo@10.163.183.71:~/Desktop/B71v2sh_masked.fasta blast/ 
cd blast
blastn -query B71v2sh_masked.fasta -subject UFVPY184_nh.fasta -evalue 1e-50 -max_target_seqs 20000 -outfmt '6 qseqid sseqid qstart qend sstart send btop' -out B71v2sh.UFVPY184.BLAST
```

## 21. Identify SNP variants between UFVPY184 and B71v2sh genome. 
```bash
ssh cjea222@mcc.uky.edu
cd /project/farman_s24cs485g/cjea222/
scp cjea222@cjea222.cs.uky.edu:~/blast/B71v2sh.UFVPY184.BLAST
cd ..
cp SCRIPTs/CallVariants.sh cjea222/
cd cjea222
sbatch CallVariants.sh B71v2sh.UFVPY184.BLAST
```

## 22. Remove any contigs in UFVPY184 less than 200 base pairs in length.
```bash
ssh cjea222@mcc.uky.edu
cd /project/farman_s24cs485g/
cp SCRIPTs/CullShortContigs.pl cjea222/
cd cjea222
perl CullShortContigs.pl UFVPY184_nh.fasta
```
