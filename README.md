# RNA-MosaicHunter

## 1. Software requirement
- java 1.8 (https://java.com/en/)  
- UCSC blat (http://genome.ucsc.edu/FAQ/FAQblat.html#blat3)  
- STAR 2.5.4
- Picard 1.138
- GATK 3.6
- SAMtools 0.1.19
- BEDTools 2.23.0


## 2. Prepare reference
RNA-MosaicHunter requires a fasta file (.fasta, .fa) for your reference genome. It is better to make sure that your reference file has the same name and order of contigs to your .bam file(s) and .bed file(s). Reads aligned to any contigs which do not appear in the reference file will be ignored. When running RNA-MosaicHunter, you need to set the top parameter reference_file.

## 3. Prepare data for RNA-MosaicHunter Input

RNA-MosaicHunter currently accepts aligned reads from RNA sequencing (RNA-seq). We recommend STAR for read mapping and GATK/Picard for read pre-processing. These pre-processing steps are important in order to reduce false positives in identifying mosaic sites. The following is an example of how the reads are pre-processed.  
  
Note that default config file and resources files are designed for GRCh37 (hg19). To run on other versions (e.g. GRCh38), users should prepare hg38 reference files and change the parameters "valid_references", "chr_x_name", "chr_y_name" in the config file to the corresponding chromosome names (e.g. with "chr" for GRCh38).  


### 3.1 Alignment -- STAR
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

### 3.4 Local realignment -- GATK
```
java -jar GenomeAnalysisTK.jar -T RealignerTargetCreator -nt ${THREAD_NUM} -R ${REFERENCE_DIR}/human_v37_contig_hg19_hs37d5/human_v37_contig_hg19_hs37d5.fasta -I rnaseq/GATK/${IND_NAME}_${SAMPLE_ID}.split.bam -o rnaseq/GATK/${IND_NAME}_${SAMPLE_ID}.intervals

java -jar GenomeAnalysisTK.jar -T IndelRealigner -R ${REFERENCE_DIR}/human_v37_contig_hg19_hs37d5/human_v37_contig_hg19_hs37d5.fasta -I rnaseq/GATK/${IND_NAME}_${SAMPLE_ID}.split.bam -targetIntervals rnaseq/GATK/${IND_NAME}_${SAMPLE_ID}.intervals -o rnaseq/GATK/${IND_NAME}_${SAMPLE_ID}.realigned.bam
echo "Job ended at $(date)"
```

### 3.5 Read calibration -- GATK
```
java -jar GenomeAnalysisTK.jar -T BaseRecalibrator -R ${REFERENCE_DIR}/human_v37_contig_hg19_hs37d5/human_v37_contig_hg19_hs37d5.fasta -I rnaseq/GATK/${IND_NAME}_${SAMPLE_ID}.realigned.bam -knownSites ${REFERENCE_DIR}/dbsnp_137.b37.vcf -o rnaseq/GATK/${IND_NAME}_${SAMPLE_ID}.recal_data.grp

java -jar GenomeAnalysisTK.jar -T PrintReads -nct ${THREAD_NUM} -R ${REFERENCE_DIR}/human_v37_contig_hg19_hs37d5/human_v37_contig_hg19_hs37d5.fasta -I rnaseq/GATK/${IND_NAME}_${SAMPLE_ID}.realigned.bam -BQSR rnaseq/GATK/${IND_NAME}_${SAMPLE_ID}.recal_data.grp -o rnaseq/GATK/${IND_NAME}_${SAMPLE_ID}.BQSR.bam
```

### 3.6 Remove improper-paired and multiple-hit reads -- Samtools
```
samtools view -f 2 -F 768 -h rnaseq/GATK/${IND_NAME}_${SAMPLE_ID}.BQSR.bam | samtools view -Sb - > rnaseq/bam_ready/${IND_NAME}_${SAMPLE_ID}.final.bam
samtools index rnaseq/bam_ready/${IND_NAME}_${SAMPLE_ID}.final.bam
```
${IND_NAME}_${SAMPLE_ID}.final.bam could be used as the input for Step 5 Run RNA-MosaicHunter.

### 3.7 Calculate coverage -- Bedtools (not necessary)
```
/BEDtools-2.23.0/bin/coverageBed -abam rnaseq/bam_ready/${IND_NAME}_${SAMPLE_ID}.final.bam -b ${REFERENCE_DIR}/human_gene_gencode19.b37.exon.bed -hist | awk '$1=="all"' > rnaseq/coverageBed/${IND_NAME}_${SAMPLE_ID}.tsv
```
  
## 4. Example data for RNA-MosaicHunter Input
resources/example.bam (generated based on sequence around GTEX-111FC (sSNV chr4:6064082)). Need preprocessing including sorting and indexing.
```
samtools sort example.bam example_srt
samtools index example_srt.bam
```

## 5. Run RNA-MosaicHunter
### 5.1 RNA-MosaicHunter
```
java -Xmx64G -jar ${MOSAICHUNTER_DIR}/build/mosaichunter.jar -C ${MOSAICHUNTER_CONFIG} -P input_file=${PATH}/example_srt.bam -P mosaic_filter.sex=${SEX} -P output_dir=${OUTPUT_DIR} -P common_site_filter.bed_file=${ERROR_PRONE_BED} -P misaligned_reads_filter.max_NM=3

cat ${OUTPUT_DIR}/final.passed.tsv | awk '$3==$7||$3==$9' | awk '$11~"N/A"&&$11~"1.0"' > ${OUTPUT_DIR}/final.clean.tsv
```
Example Output of RNA-MosaicHunter:
```
base_quality = 33
chr_x_name = X
chr_y_name = Y
depth_sampling = true
enable_reference_cache = false
in_process_filter_name = in_process_filter
input_file = test_0122s_srt.bam
input_sampling = false
input_sampling_regions = 1
input_sampling_size = 1
max_depth = 5001
max_recent_reads = 9997
min_mapping_quality = 20
min_read_quality = 20
output_dir = tmp_0122
post_process_filter_name = post_process_filter
read_buffer_size = 100000
reference_file = /home/yh174/reference/human_v37/human_g1k_v37.fasta
remove_duplicates = true
remove_flags = 0x100
seed = 0
valid_references = 1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,X,Y
base_number_filter.class = cn.edu.pku.cbi.mosaichunter.filter.BaseNumberFilter
base_number_filter.min_minor_allele_number = 5
base_number_filter.min_minor_allele_percentage = 5
base_number_filter.omit_alt_homozygous = false
base_number_filter.output_passed = true
common_site_filter.bed_file = /home/yh174/reference/error_prone_bed/FSJ-M_XHX-F_ACC1-blood_PE.overlap.singularAF.bed
complete_linkage_filter.binom_error_rate = 1e-3
complete_linkage_filter.binom_p_value_cutoff = 0.01
complete_linkage_filter.class = cn.edu.pku.cbi.mosaichunter.filter.CompleteLinkageFilter
complete_linkage_filter.fisher_p_value_cutoff = 0.01
complete_linkage_filter.output_filtered = true
complete_linkage_filter.output_passed = true
depth_filter.class = cn.edu.pku.cbi.mosaichunter.filter.DepthFilter
depth_filter.max_depth = 500000
depth_filter.min_depth = 10
depth_filter.output_filtered = false
depth_filter.output_passed = false
final.class = cn.edu.pku.cbi.mosaichunter.filter.OutputFilter
final.data_name = mosaic_filter
final.output_passed = true
homopolymers_filter.class = cn.edu.pku.cbi.mosaichunter.filter.HomopolymersFilter
homopolymers_filter.long_homopolymer_expansion = 3
homopolymers_filter.long_homopolymer_length = 6
homopolymers_filter.output_filtered = true
homopolymers_filter.output_passed = true
homopolymers_filter.short_homopolymer_expansion = 2
homopolymers_filter.short_homopolymer_length = 4
in_process_filter.class = cn.edu.pku.cbi.mosaichunter.filter.AndFilter
in_process_filter.filters = depth_filter,base_number_filter,repetitive_region_filter,homopolymers_filter,strand_bias_filter,mapping_quality_filter,within_read_position_filter,mosaic_filter,complete_linkage_filter
mapping_quality_filter.class = cn.edu.pku.cbi.mosaichunter.filter.MappingQualityFilter
mapping_quality_filter.output_filtered = true
mapping_quality_filter.output_passed = true
mapping_quality_filter.p_value_cutoff = 0.05
misaligned_reads_filter.blat_param = -stepSize=5 -repMatch=2253 -minScore=0 -minIdentity=0.5 -noHead
misaligned_reads_filter.blat_path = 
misaligned_reads_filter.class = cn.edu.pku.cbi.mosaichunter.filter.MisalignedReadsFilter
misaligned_reads_filter.max_NM = 3
misaligned_reads_filter.min_gap_distance = 5
misaligned_reads_filter.min_overlap_percentage = 0.9
misaligned_reads_filter.min_side_distance = 5
misaligned_reads_filter.misalignment_threshold = 0.5
misaligned_reads_filter.omit_multiple_alignment = false
misaligned_reads_filter.output_filtered = true
misaligned_reads_filter.output_passed = true
misaligned_reads_filter.reference_file = /home/yh174/reference/human_v37_contig_hg19_hs37d5.fasta
mosaic_filter.alpha_param = 0
mosaic_filter.base_change_rate = 1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1
mosaic_filter.beta_param = 0
mosaic_filter.class = cn.edu.pku.cbi.mosaichunter.filter.MosaicFilter
mosaic_filter.control_bam_file = 
mosaic_filter.dbsnp_file = /home/yh174/reference/dbsnp_137.b37.raw.tsv
mosaic_filter.de_novo_rate = 1e-8
mosaic_filter.father_bam_file = 
mosaic_filter.fisher_threshold = 0.01                   
mosaic_filter.mode = RNA
mosaic_filter.mosaic_rate = 1e-7
mosaic_filter.mosaic_threshold = 0.05
mosaic_filter.mother_bam_file = 
mosaic_filter.novel_af = 1e-4
mosaic_filter.omit_alt_homozygous = false
mosaic_filter.output_filtered = true
mosaic_filter.output_passed = true
mosaic_filter.sex = M
mosaic_filter.unknown_af = 0.002
null_filter.class = cn.edu.pku.cbi.mosaichunter.filter.NullFilter
null_filter.output_filtered = true
null_filter.output_passed = true
null_filter.return_value = false
post_process_filter.class = cn.edu.pku.cbi.mosaichunter.filter.AndFilter
post_process_filter.filters = misaligned_reads_filter,final
repetitive_region_filter.bed_file = /home/yh174/reference/repeats/all_repeats.b37.bed
repetitive_region_filter.class = cn.edu.pku.cbi.mosaichunter.filter.RegionFilter
repetitive_region_filter.include = false
repetitive_region_filter.output_filtered = true
repetitive_region_filter.output_passed = true
stats_manager.enable_counter = false
stats_manager.enable_timer = false
strand_bias_filter.class = cn.edu.pku.cbi.mosaichunter.filter.StrandBiasFilter
strand_bias_filter.output_filtered = true
strand_bias_filter.output_passed = true
strand_bias_filter.p_value_cutoff = 0.05
within_read_position_filter.class = cn.edu.pku.cbi.mosaichunter.filter.WithinReadPositionFilter
within_read_position_filter.output_filtered = true
within_read_position_filter.output_passed = true
within_read_position_filter.p_value_cutoff = 0.05
Wed Jan 22 01:19:49 EST 2025 Initializing...
Wed Jan 22 01:19:51 EST 2025 Reading reference from file: /home/yh174/reference/human_v37/human_g1k_v37.fasta
Wed Jan 22 01:20:25 EST 2025 Initializing filters...
Wed Jan 22 01:20:29 EST 2025 Scanning...
Wed Jan 22 01:20:29 EST 2025 - Time(s):0 Reads:0 Sites:0/3095677412 Progress:0.00%
Wed Jan 22 01:20:40 EST 2025 - Time(s):10 Reads:127 Sites:3095677412/3095677412 Progress:100.00%
run blat: blat -stepSize=5 -repMatch=2253 -minScore=0 -minIdentity=0.5 -noHead /home/yh174/reference/human_v37_contig_hg19_hs37d5.fasta tmp_0122/misaligned_reads_filter.fa tmp_0122/misaligned_reads_filter.psl
filter name                                          pass/all   ratio
depth_filter                                          111/144  77.08%
base_number_filter                                      1/111   0.90%
repetitive_region_filter                                  1/1 100.00%
homopolymers_filter                                       1/1 100.00%
strand_bias_filter                                        1/1 100.00%
mapping_quality_filter                                    1/1 100.00%
within_read_position_filter                               1/1 100.00%
mosaic_filter                                             1/1 100.00%
complete_linkage_filter                                   1/1 100.00%
misaligned_reads_filter                                   1/1 100.00%
final                                                     1/1 100.00%
```


### 5.2 RNA-editing filter
We have two different mode of RNA editint filters:
1. **full-removal**: removes all A>G / T>C mutations  

```
Rscript RNA_Editing_Filter_MH.R ${OUTPUT_DIR}/final.clean.tsv full-removal ${OUTPUT_DIR}/final.clean.filter_RNA_edit.tsv
```

2. **filter-based-removal**: removes RNA editing sites reported in DARNED [1] and REDIportal [2], and A>G mutations on transcribed strand and T>C mutations on untranscribed strand  

```
Rscript RNA_Editing_Filter_MH.R ${OUTPUT_DIR}/final.clean.tsv filter-based-removal ${OUTPUT_DIR}/final.clean.filter_RNA_edit.tsv
```

## 6. Output format

All of the output files are in the specified directory output_dir, including the final candidate list final.passed.tsv, the filtered and remained remaining variant lists of each filters (.filtered.tsv, .passed.tsv), and other temporary files, such as blat input and output.


The final.clean.filter_RNA_edit.tsv files are in tab-separated format, whose columns' meanings for single mode are listed below:  

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

## References

[1] Kiran,A. and Baranov,P.V. (2010) DARNED: a DAtabase of RNa EDiting in humans. Bioinformatics, 26, 1772–1776.  
[2] Picardi,E., D’Erchia,A.M., Lo Giudice,C. and Pesole,G. (2017) REDIportal: a comprehensive database of A-to-I RNA editing events in humans. Nucleic Acids Res., 45, D750–D757.  

