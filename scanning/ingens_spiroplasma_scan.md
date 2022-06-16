# Scanning *Spiroplasma* from the samples of *N. ingens*
The script is used for scanning genes of the *Spiroplasma* from the low-coverage whole genome sequences data of the *N. ingens* samples.

1. map the raw reads of the beetles to the *Spiroplasma* genome
- create a list to run bwa (a mapping approach)
``` 
for file in /media/Jade/ingens_low-coverage_seq/raw_fq/*/*.gz
do
    echo ${file} | cut -d "/" -f 7 | cut -d "_" -f 1,2 >> list
done
cat list | sort | uniq >> list2
rm list
mv list2 list
```

- run bwa with the created list
```
bwa index ./spiro_ref_genome/spiroplasma.txt

# divide the list into list1, list2...list6 to paralyze the process with multiple threads.
# here only the example of list 6 is shown
for file in $(cat ./list6)
do
    echo "working on ${file}...."
    R1=$(ls /media/Jade/ingens_low-coverage_seq/raw_fq/*/${file}* | grep "R1")
    R2=$(ls /media/Jade/ingens_low-coverage_seq/raw_fq/*/${file}* | grep "R2")
    bwa mem ./spiro_ref_genome/spiroplasma.txt ${R1} ${R2} > ${file}.sam
    samtools view -S -b ${file}.sam > ${file}.bam
    rm ${file}.sam
done
```

2. Calculate the coverage of each bam files
```
# working dir: sean@Fuji:/media/Jade/YMW/spiroscan/
# create a list to make the unit in the base of individual beetle
# This is done by merging bam files from same individual

ls *.bam | cut -d "_" -f 1 | sort | uniq >> beetle_list
for beetle in $(cat beetle_list)
do
    count=$(ls *.bam | grep "${beetle}" | cut -d "_" -f 1 | grep -ow "${beetle}" | wc -l)
    echo -e "${beetle}\t${count}" >> beetle_with_count
done
cat beetle_with_count | grep -v "1$" >> duplicated_sample 
# there are 191 samples with at least 2 bam files

# get the mapped read and write into new bam file
mkdir mapped_reads
for file in *bam
do
    ../../../Data1/Yiming/software/samtools-1.13/samtools view -b -F 4 ${file} > ./mapped_reads/${file}
done

# sort the mapped bam files
# working dir: sean@Fuji:/media/Jade/YMW/spiroscan/mapped_reads
for file in *.bam
do
    name=$(echo ${file} | cut -d "." -f 1)
    ../../../../Data1/Yiming/software/samtools-1.13/samtools sort ${file} -o ${name}_sorted.bam
done

# merge bam files according to the duplicate list, with samtools
IFS=$'\n'
mkdir tmp
for beetle in $(cat ../duplicated_sample | cut -d $'\t' -f 1)
do
    echo "working on ${beetle}..."
    ../../../../Data1/Yiming/software/samtools-1.13/samtools merge -o ${beetle}_merged_sorted.bam ${beetle}_*.bam
    mv ${beetle}_S*.bam tmp
done

### remove original bwa output bam files, only keep the mapped_sorted bam files for calculating the coverage and depth

# calculate the coverage and depth for the bam files
mkdir coverage
for file in *.bam
do
   name=$(echo $file | cut -d "_" -f 1)
   ../../../../Data1/Yiming/software/bbmap/pileup.sh in=${file}  out=./coverage/${name}.out
   covered=$(cat ./coverage/${name}.out | grep "tig00009077" | cut -d $'\t' -f 5)
   echo -e "${name}\t${covered}" >> all_coverage
done
```
  - defining the infection is the sample with coverage >= 70%, the overall infection rate is 56.14% 
  ![](@attachment/Clipboard_2021-08-03-00-32-01.png)
  - These sample were not considered due to the low sequencing depth:
  **SEKI-0747E, YMW18-019, YMW18-025, YMW18-031, YMW18-041, YMW18-043, YMW18-048, YMW18-071**
  - <span style="color:blue"> the original bam output from bwa have been removed </span>
  - <span style="color:blue"> the mapped but unsorted bam files also have been removed </span>.

3. visualize the distribution of genome coverage rate (the percentage of the *Spiroplasma* genome been covered by from the mapped reads) with R.
- The R script is written in an independent script file called "scan_plot.r".

```
setwd("C:/Users/wengz/Box/Spiroplasma/scanning")

coverage <- read.table("all_coverage", header=F, sep="\t")
hist(coverage$V2, main = "covered_percent", xlab= "percentage", ylab="count")
```