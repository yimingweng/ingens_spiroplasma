# 2021.08.02-calling variants from Spiroplasma bam files
##  This is for making the population structure of the spiroplasma
1. Calling SNPs with GATK
2. To call the variats, I follow the [GATK pipeline](https://gatk.broadinstitute.org/hc/en-us/articles/360035535932-Germline-short-variant-discovery-SNPs-Indels-)
3. With the sorted bam file, I will need to mark the duplicates using [MarkDuplicates (Picard)](https://gatk.broadinstitute.org/hc/en-us/articles/360037052812-MarkDuplicates-Picard-):
```
# sean@Fuji:/media/Jade/YMW/spiroscan/vcf
mkdir MarkDuplicates
# Use Picard to mark the duplicates for each sort-mapped bam files
for file in ../mapped_reads/infected/*.bam
do
  name=$(echo ${file} | cut -d "/" -f 4 | cut -d "_" -f 1)
  echo "working ${name}..."
  java -jar /media/Data1/Yiming/software/picard.jar MarkDuplicates I=${file} O=./MarkDuplicates/${name}_markdups.bam M=./MarkDuplicates/${name}_markdups.txt
done
```

4. Because the standard variants are not available for this data to do BQSR, I will need to generate a variants calling data without Base Quality Score Recalibration [(BQSR)](https://gatk.broadinstitute.org/hc/en-us/articles/360035890531-Base-Quality-Score-Recalibration-BQSR-), and use that VCF as standard to run BQSR.
5. So do the haplotypecaller in GVCF form first:
```
# fai index the genome for the reference genome
# working dir: sean@Fuji:/media/Jade/YMW/spiroscan/spiro_ref_genome
# here I changed the name of the genome from spiroplasma.txt to spiroplasma.fasta
/media/Data1/Yiming/software/samtools-1.13/samtools faidx spiroplasma.fasta

# add @RG for the bam file (require for haplotypecaller)
# working dir: sean@Fuji:/media/Jade/YMW/spiroscan/vcf/MarkDuplicates/
mkdir ../RGbam
for file in *.bam
do
  name=$(echo "${file}"| cut -d "_" -f 1)
  echo "working on ${name}..."
  java -jar /media/Data1/Yiming/software/picard.jar AddOrReplaceReadGroups INPUT=${file} OUTPUT=../RGbam/${name}.RG.bam RGID=ingens RGLB=library1 RGPL=illumina RGPU=unit1 RGSM=${name}
done

# working dir: sean@Fuji:/media/Jade/YMW/spiroscan/vcf/RGbam/
# create sequence dir for the reference genome (not sure this is required)
/media/Data1/Yiming/software/gatk-4.1.9.0/gatk CreateSequenceDictionary -R /media/Jade/YMW/spiroscan/spiro_ref_genome/spiroplasma.fasta

# index all .RG.bam files
for file in *bam
do
  /media/Data1/Yiming/software/samtools-1.13/samtools index ${file}
done

# run haplotypecaller and write output into ./GVCF/
# To parallize the process, make 7 lists and each takes 30 samples of the 210.
# Here I only show the example for list1
mkdir ../GVCF
ls *.bam | head -30 >> list1
ls *.bam | head -60 | tail -30 >> list2
IFS=$'\n'
for file in $(cat list7)
do
   name=$(echo $file | cut -d "." -f 1)
   /media/Data1/Yiming/software/gatk-4.1.9.0/gatk --java-options "-Xmx15g" HaplotypeCaller \
   -R ../../spiro_ref_genome/spiroplasma.fasta \
   -I  ${file} \
   -O ../GVCF/${name}.g.vcf.gz \
   -ploidy 1 \
   -ERC GVCF
done
```
2. Combind GVCF files from haplotypecaller into single VCF file with [GenomicsDBImport](https://gatk.broadinstitute.org/hc/en-us/articles/360036883491-GenomicsDBImport)
```
# working dir: sean@Fuji:/media/Jade/YMW/spiroscan/vcf/
# create output dir
mkdir GenomicsDBImport

# create a map file contain all the inpute GVCF files
# working in /data2/ingens_low_coverage/GVCF/
for file in *.vcf.gz
do 
  name=$(echo $file | cut -d "." -f 1) 
  path=$(echo "/media/Jade/YMW/spiroscan/vcf/GVCF/${file}")
  echo -e "${name}\t${path}" >> ../GenomicsDBImport/gvcfmap
done

# create contig list, even there are only 2 contigs
cat /media/Jade/YMW/spiroscan/spiro_ref_genome/spiroplasma.fasta | grep ">" | cut -d " " -f 1 | sed 's/>//g' >> /media/Jade/YMW/spiroscan/vcf/GenomicsDBImport/contiglist

# run GenomicsDBImport for each contig
mkdir tmp
for contig in $(cat contiglist)
do
  /media/Data1/Yiming/software/gatk-4.1.9.0/gatk --java-options "-Xmx32g -Xms32g" GenomicsDBImport --genomicsdb-workspace-path DBImport_${contig} -L ${contig} --sample-name-map gvcfmap --tmp-dir ./tmp --reader-threads 28
done

# run GenotypeGVCFs for each contig
for contig in $(cat contiglist)
do
  /media/Data1/Yiming/software/gatk-4.1.9.0/gatk --java-options "-Xmx32g" GenotypeGVCFs -R /media/Jade/YMW/spiroscan/spiro_ref_genome/spiroplasma.fasta -V gendb://DBImport_${contig} -O ${contig}.vcf.gz
done
```