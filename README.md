# RNA-MosaicHunter

## 1. Software requirement
- java 1.8 (https://java.com/en/)  
- UCSC blat (http://genome.ucsc.edu/FAQ/FAQblat.html#blat3)  
- STAR 2.5.4
- Picard 1.138
- GATK 3.6
- SAMtools 0.1.19
- BEDTools 2.26.0

## 2. Prepare reference
RNA-MosaicHunter requires a `fasta` file (.fasta, .fa) for your reference genome. It is better to make sure that your reference file has the same name and order of contigs to your `.bam` file(s) and `.bed` file(s). Reads aligned to any contigs that do not appear in the reference file will be ignored. When running RNA-MosaicHunter, you need to set the top parameter reference_file.

## 3. Prepare data for RNA-MosaicHunter input

RNA-MosaicHunter currently accepts aligned reads from RNA sequencing (RNA-seq). We recommend STAR for read mapping and GATK/Picard for read pre-processing. These pre-processing steps are important in order to reduce false positives in identifying mosaic sites. The following is an example of how the reads are pre-processed.  
  
Note that default config file and resources files are designed for GRCh37 (hg19). To run on other versions (e.g. GRCh38), users should prepare hg38 reference files and change the parameters `valid_references`, `chr_x_name`, `chr_y_name` in the config file to the corresponding chromosome names (e.g. with "chr" for GRCh38).  

### 3.1 Read alignment -- STAR
```
STAR --runThreadN ${THREAD_NUM} --twopassMode Basic --outSAMattributes All --outSAMtype BAM Unsorted --readFilesCommand zcat --genomeDir ${REFERENCE_DIR}/human_v37_contig_hg19_hs37d5_STAR_index10_gencode19 --readFilesIn "rnaseq/fastq/${IND_NAME}/${SAMPLE_ID}_1.fastq.gz" "rnaseq/fastq/${IND_NAME}/${SAMPLE_ID}_2.fastq.gz" --outFileNamePrefix "$OUTPUT_DIR/${IND_NAME}/${SAMPLE_ID}/"
```

### 3.2 Remove duplicates -- Picard
```
java -jar picard.jar AddOrReplaceReadGroups INPUT=rnaseq/STAR/${IND_NAME}/${SAMPLE_ID}/Aligned.out.bam OUTPUT=rnaseq/Picard/${IND_NAME}_${SAMPLE_ID}.sorted.bam SO=coordinate ID=${IND_NAME}_${SAMPLE_ID} LB=unknown PL=illumina SM=${IND_NAME}_${SAMPLE_ID} PU=unknown
java -jar picard.jar MarkDuplicates INPUT=rnaseq/Picard/${IND_NAME}_${SAMPLE_ID}.sorted.bam OUTPUT=rnaseq/Picard/${IND_NAME}_${SAMPLE_ID}.masked.bam M=rnaseq/Picard/${IND_NAME}_${SAMPLE_ID}.matrix REMOVE_DUPLICATES=true AS=true VALIDATION_STRINGENCY=SILENT
samtools index rnaseq/Picard/${IND_NAME}_${SAMPLE_ID}.masked.bam
```

### 3.3 Split'N'Trim and mapping quality reassignment -- GATK
```
java -jar GenomeAnalysisTK.jar -T SplitNCigarReads -R ${REFERENCE_DIR}/human_v37_contig_hg19_hs37d5/human_v37_contig_hg19_hs37d5.fasta -I rnaseq/Picard/${IND_NAME}_${SAMPLE_ID}.masked.f.bam -o rnaseq/GATK/${IND_NAME}_${SAMPLE_ID}.split.bam -rf ReassignOneMappingQuality -RMQF 255 -RMQT 60 -U ALLOW_N_CIGAR_READS
```

### 3.4 Local indel realignment -- GATK
```
java -jar GenomeAnalysisTK.jar -T RealignerTargetCreator -nt ${THREAD_NUM} -R ${REFERENCE_DIR}/human_v37_contig_hg19_hs37d5/human_v37_contig_hg19_hs37d5.fasta -I rnaseq/GATK/${IND_NAME}_${SAMPLE_ID}.split.bam -o rnaseq/GATK/${IND_NAME}_${SAMPLE_ID}.intervals
java -jar GenomeAnalysisTK.jar -T IndelRealigner -R ${REFERENCE_DIR}/human_v37_contig_hg19_hs37d5/human_v37_contig_hg19_hs37d5.fasta -I rnaseq/GATK/${IND_NAME}_${SAMPLE_ID}.split.bam -targetIntervals rnaseq/GATK/${IND_NAME}_${SAMPLE_ID}.intervals -o rnaseq/GATK/${IND_NAME}_${SAMPLE_ID}.realigned.bam
```

### 3.5 Base quality calibration -- GATK
```
java -jar GenomeAnalysisTK.jar -T BaseRecalibrator -R ${REFERENCE_DIR}/human_v37_contig_hg19_hs37d5/human_v37_contig_hg19_hs37d5.fasta -I rnaseq/GATK/${IND_NAME}_${SAMPLE_ID}.realigned.bam -knownSites ${REFERENCE_DIR}/dbsnp_137.b37.vcf -o rnaseq/GATK/${IND_NAME}_${SAMPLE_ID}.recal_data.grp
java -jar GenomeAnalysisTK.jar -T PrintReads -nct ${THREAD_NUM} -R ${REFERENCE_DIR}/human_v37_contig_hg19_hs37d5/human_v37_contig_hg19_hs37d5.fasta -I rnaseq/GATK/${IND_NAME}_${SAMPLE_ID}.realigned.bam -BQSR rnaseq/GATK/${IND_NAME}_${SAMPLE_ID}.recal_data.grp -o rnaseq/GATK/${IND_NAME}_${SAMPLE_ID}.BQSR.bam
```

### 3.6 Remove improper-paired and multiple-hit reads -- Samtools
```
samtools view -f 2 -F 768 -h rnaseq/GATK/${IND_NAME}_${SAMPLE_ID}.BQSR.bam | samtools view -Sb - > rnaseq/bam_ready/${IND_NAME}_${SAMPLE_ID}.final.bam
samtools index rnaseq/bam_ready/${IND_NAME}_${SAMPLE_ID}.final.bam
```
`${IND_NAME}_${SAMPLE_ID}.final.bam` could be used as the input for Step 4 Run RNA-MosaicHunter.

## 4. Run RNA-MosaicHunter
### 4.1 RNA-MosaicHunter
```
java -Xmx64G -jar ${MOSAICHUNTER_DIR}/build/mosaichunter.jar -C ${MOSAICHUNTER_CONFIG} -P input_file=${PATH}/example_srt.bam -P mosaic_filter.sex=${SEX} -P output_dir=${OUTPUT_DIR}
```

You can modify RNA-MosaicHunter parameters either by editing the `${MOSAICHUNTER_CONFIG}` file or by adding `-P ${PARAMETER_NAME}=${PARAMETER_VALUE}` to the command line.
Default version of `${MOSAICHUNTER_CONFIG}` based on GRCh37 is provided under `/conf`.

### 4.2 Remove sites where the top two alleles are both non-reference or those present in dbSNP (optional)
```
cat ${OUTPUT_DIR}/final.passed.tsv | awk '$3==$7||$3==$9' | awk '$11~"N/A"&&$11~"1.0"' > ${OUTPUT_DIR}/final.clean.tsv
```

### 4.3 RNA-editing filter
We have two different mode of RNA editint filters:
1. **full-removal**: removes all A>G / T>C mutations  

```
Rscript RNA_Editing_Filter_MH.R ${OUTPUT_DIR}/final.clean.tsv full-removal ${OUTPUT_DIR}/final.clean.filter_RNA_edit.tsv
```

2. **filter-based-removal**: removes RNA editing sites reported in DARNED [1] and REDIportal [2], and A>G mutations on transcribed strand and T>C mutations on untranscribed strand  

```
Rscript RNA_Editing_Filter_MH.R ${OUTPUT_DIR}/final.clean.tsv filter-based-removal ${OUTPUT_DIR}/final.clean.filter_RNA_edit.tsv
```

## 5. Output format

All of the output files are in the specified directory ${OUTPUT_DIR}, including the final mosaic candidate list `final.clean.filter_RNA_edit.tsv`, the filtered and remained remaining variant lists of each filters (`.filtered.tsv`, `.passed.tsv`), and other temporary files, such as blat input and output.


The final.clean.filter_RNA_edit.tsv files are in tab-separated format, with column meanings detailed below:  

1. Contig / chromosome name
2. Position / coordinate on the contig (1-based)
3. Base of reference allele
4. Total depth of this site
5. Pileuped sequencing bases at this site
6. Pileuped sequencing baseQs at this site
7. Base of major allele
8. Depth of major allele
9. Base of minor allele
10. Depth of minor allele
11. dbSNP allele frequency of major and minor alleles
12. log10 prior probability of major-homozygous genotype
13. log10 prior probability of heterozygous genotype
14. log10 prior probability of minor-homozygous genotype
15. log10 prior probability of mosaic genotype
16. log10 likelihood of major-homozygous genotype
17. log10 likelihood of heterozygous genotype
18. log10 likelihood of minor-homozygous genotype
19. log10 likelihood of mosaic genotype
20. log10 posterior probability of major-homozygous genotype
21. log10 posterior probability of heterozygous genotype
22. log10 posterior probability of minor-homozygous genotype
23. log10 posterior probability of mosaic genotype
24. Mosaic posterior probability

You can get a high-confidence candidate list of mosaic sites by filtering the 24th column (mosaic posterior probability) or sorting this column from high to low.
  
## 6. Demo
The example input of RNA-MosaicHunter can be found at `resources/example.bam`:
```
samtools sort example.bam example_srt
samtools index example_srt.bam
java -Xmx64G -jar build/mosaichunter.jar -C conf/RNA.properties -P input_file=example_srt.bam -P mosaic_filter.sex=F -P output_dir=. -P misaligned_reads_filter.max_NM=3
cat final.passed.tsv | awk '$3==$7||$3==$9' | awk '$11~"N/A"&&$11~"1.0"' > final.clean.tsv
Rscript RNA_Editing_Filter_MH.R final.clean.tsv filter-based-removal final.clean.filter_RNA_edit.tsv
```

The expected output of RNA-MosaicHunter (`final.clean.filter_RNA_edit.tsv`):
```
4	6064082	TRUE	122	TATTTTTTTTTTTTTTATTTTttttttttttttttttttttttttttttttttTTtttttTTTAttTttttatttttTTTTTAtTttTTTTTTAtTtatttttTtttttatatttttTttTt	BDCACBCBBBBCBBBBCBBBBCBCCCCCCCBCCCCDCCCDCCCBCCACCCCACBBCCCCCCBBCCCBCCC>BCCA>CACBCCDDBCDBBBBBBCDADBDDDDDBECDDDBABDDDDD=DDCH	TRUE	113	A	9	T(1.0):A(N/A-9e-05	-3.69906	-8.00009	-7	-30.42197	-36.72566	-380.90166	-15.17676	-8.2453	-18.24796	-366.72499	0	1
```

## References

[1] Kiran A and Baranov PV (2010) DARNED: a DAtabase of RNa EDiting in humans. Bioinformatics, 26:1772–1776.  
[2] Picardi E, D’Erchia AM, Lo Giudice C, Pesole G (2017) REDIportal: a comprehensive database of A-to-I RNA editing events in humans. Nucleic Acids Res, 45:D750–D757.  

