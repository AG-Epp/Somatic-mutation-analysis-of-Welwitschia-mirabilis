# Somatic mutation analysis of Welwitschia mirabilis
This repository contains digital supplements for the master thesis of Clara T., "Somatic mutation analysis for age estimation of Welwitschia mirabilis".

The data processing for this project includes

	- mapping of Illumina sequencing data from three Welwitschia individuals using paleomix (2 samples per individual)
 
 	- variant calling using GATK HaplotypeCaller
  
  	- counting intra-individual somatic variants using a custom awk script
   
   	- analysis of variant counts and inter-variant distances using R
    
It is currently not possible to include the entirety of the files and data to run this analysis in full, however intermediate files are given for R processing, and the main data processing and file usage is detailed below.

# Bioinformatic processing for joint variant calling using GATK

Paleomix mapping

	# yaml file (include file)
	
	Welwitschia_mapping.2024-03-01.yaml	
	
	# paleomix command line
	
		paleomix bam_pipeline run /.../Welwitschia_mapping.2024-03-01.yaml 
		--jar-root=/data/vanessa/jar_root 
		--max-threads 100 
		--bwa-max-threads 12 
		--adapterremoval-max-threads 24
		--temp-root /data/lastexpansion/clarata/temp/

# BAM to VCF file workflow
Command lines are exemplary.
SAMPLE is placeholder for samples, i.e. WW01_1, WW01_2, WW04_1, WW04_5, ME01_1 or ME01_3
CHR is placeholder for chromosome, i.e. Chr01, Chr02, ..., Chr21
commands were executed for every chromosome of every sample.


# BAM file sorting 
	# Create single chromosome BAM files for every sample
	
	samtools view -b SAMPLE.welwitschia_mirabilis_nuclear.bam CHR > SAMPLE.welwitschia_mirabilis_nuclear.CHR.bam
		
	# Sort BAM files
	
	samtools sort SAMPLE.welwitschia_mirabilis_nuclear.CHR.bam > SAMPLE.welwitschia_mirabilis_nuclear.CHR.sort.bam 
	
	# Index BAM file 
	
	samtools index -@ 4 -c SAMPLE.welwitschia_mirabilis_nuclear.CHR.sort.bam 
	
	
# Variant calling with GATK HaplotypeCaller, creating genomic VCF files
	
	java -jar /.../gatk-package-4.3.0.0-local.jar HaplotypeCaller 
		-R /.../Welwitschia_genome.fasta 
		--min-base-quality-score 20 
		--minimum-mapping-quality 30 
		-ERC BP_RESOLUTION 
		-ploidy 2 
		--output-mode EMIT_ALL_ACTIVE_SITES 
		-I /.../SAMPLE.welwitschia_mirabilis_nuclear.CHR.sort.bam 
		-O SAMPLE.welwitschia_mirabilis_nuclear.CHR.sort.bam.gatk4300.minQ20.minMAPQ30.ERC-BP_RESOLUTION.CHR.gvcf 
		--intervals CHR &
		
# Data processing for joint variant calling

	# 1. Get list of gvcf files for every chromosome
	
	ls *.CHR.gvcf > ./240306.GATK.all_samples.welwitschia.CHR_gvcf.list
		
	# 2. Run GenomicsDBimport per chromosome
	
	java -jar /.../gatk-package-4.3.0.0-local.jar GenomicsDBImport 
			-V 240306.GATK.all_samples.welwitschia.CHR_gvcf.list 
			--intervals CHR 
			--genomicsdb-workspace-path 240306.GATK.all_samples.welwitschia.CHR.genomicsDB 
			
	# 3. Joint variant calling with GenotypeGVCFs per chromosome
		This section of data processing had to be performed in multiple sections for every chromosome, since the command is very memor-intensive. 
		Files had to be merged again after.
	
	java -jar /.../gatk-package-4.3.0.0-local.jar GenotypeGVCFs 
		-R /../Welwitschia_genome.fasta
		-V gendb://240306.GATK.all_samples.welwitschia.CHR.genomicsDB 
		-stand-call-conf 0 
		-all-sites 
		-O 240306.GATK.all_samples.welwitschia.CHR.genotypeGVCFs.conf_0_allsites.vcf 
		--intervals CHR XX:XX
		
	# 4. MergeVCFs
	lists were created manually according to the sections the files were divided in
		
	java -jar /.../picard.jar MergeVcfs 
		I=CHR_GenotypeGVCFs.list 
		O=240415.GATK.all_samples.welwitschia.CHR.genotypeGVCFs.conf_0_allsites.MergeVcfs.vcf 
		
# Bedtools intersect was used to create the no_repeats dataset
	
	bedtools intersect 
		-header 
		-v 
		-a 240415.GATK.all_samples.welwitschia.CHR.genotypeGVCFs.conf_0_allsites.MergeVcfs.vcf 
		-b /.../Welwitschia_repeat_anno.Chr01.bed  
		> 240415.GATK.all_samples.welwitschia.CHR.genotypeGVCFs.conf_0_allsites.MergeVcfs.no_repeats.vcf 
		
# Using the custom AWK pairwise variant counting script
The script is specific to the samples used in the present study.
Filtering parameters are set from inside the awk script, meaning that multiple versions of the script were used for the different thresholds of filters.
The script encoded here includes the filter settings for the low-high filter thresholds. The other script versions are provided.
Output file name prefixes are specified via a variable from the command line.
	
	
Command line for script initiation:

	cat 240415.GATK.all_samples.welwitschia.CHR.genotypeGVCFs.conf_0_allsites.MergeVcfs.vcf 
		# concatenates the VCF file to pipe into AWK
	| grep -v '^#' 
		# removes header lines from VCF file
	| awk -F '\t' -v output_file=240508.RunID.variant_counting_v18_yes_bed_minRGQ40_minGQ90.CHR.allsites 
		# initiates awk, indicate file separator, define RunID, specify script version and output file name prefix
	-f variant_counting_v18_yes_bed_minRGQ40_minGQ90.txt 
		# provides AWK script file
	

The awk script produces separate variant count files for each chromosome, which are concatenated using standard unix commands before being imported into R for SMR calculation and data visualization.

In addition, the script writes bedfiles containing the position of every called intraindividual variant.


# Exploratory data analysis: extracting VCF annotation data using VariantsToTable

	java -jar /.../gatk-package-4.3.0.0-local.jar VariantsToTable 
		-V /.../240415.GATK.all_samples.welwitschia.Chr21.genotypeGVCFs.conf_0_allsites.MergeVcfs.vcf 
		-GF GQ 
		-O 240428.GATK.VarToTab.all_samples.welwitschia.Chr21.genotypeGVCFs.allsites.GF_GQ.table 
	
	java -jar /.../gatk-package-4.3.0.0-local.jar VariantsToTable 
		-V /.../240415.GATK.all_samples.welwitschia.Chr21.genotypeGVCFs.conf_0_allsites.MergeVcfs.vcf 
		-GF RGQ 
		-O 240428.GATK.VarToTab.all_samples.welwitschia.Chr21.genotypeGVCFs.allsites.GF_RGQ.table 
	
	java -jar /.../gatk-package-4.3.0.0-local.jar VariantsToTable 
		-V /.../240415.GATK.all_samples.welwitschia.Chr21.genotypeGVCFs.conf_0_allsites.MergeVcfs.vcf 
		-GF DP 
		-O 240428.GATK.VarToTab.all_samples.welwitschia.Chr21.genotypeGVCFs.allsites.GF_DP.table 
	
# final notes
Calculation of SMR, calculation of between-variant distance bins and data visualisation 
was performed in R version 4.2.2. Relevant code and files can be found in the R_scripts folder.
