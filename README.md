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
