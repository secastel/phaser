# phASER-POP
Addition to the phASER package that allows for the easy measurement of regulatory variant effects using POPulation-scale phased ASE data.

The base phASER package produces gene-level haplotypic expression per individual, phASER-POP allow gene-level haplotypic expression measurement files to be combined across individuals to produce a single haplotypic expression matrix, where each row is a gene and each column is an individual. This matrix can then be used to retrieve ASE data and calculate relevant statistics for a given set of eQTLs or any other variants of interest. See the image below for a workflow description.

![alt tag](https://raw.github.com/secastel/phaser/master/docs/phaser_pop_workflow.png)

Developed by [Stephane E. Castel](mailto:stephanecastel@gmail.com) in the [Lappalainen Lab](http://tllab.org) at the New York Genome Center and Columbia University Department of Systems Biology.

Runs on Python 2.7.x and has the following dependencies: [pandas](http://pandas.pydata.org), [Pysam](https://github.com/pysam-developers/pysam), [SciPy](http://www.scipy.org), [NumPy](http://www.numpy.org), [samtools](http://www.htslib.org), [tabix](http://www.htslib.org/doc/tabix.html)

# phaser_expr_matrix.py
## Usage
Aggregates gene-level haplotypic expression measurement files across samples to produce a single haplotypic expression matrix, where each row is a gene and each column is a sample.

## Arguments
### Required
* **--gene_ae_dir** - Directory containing all the sample-level haplotypic counts to be aggregated. Assumes each ".txt" file is an output from phaser_gene_ae.py
* **--features** - File in BED format (0 BASED COORDINATES - chr,start,stop,name) containing the features that haplotypic counts were produced for.
* **--o** - Output file.

### Optional
* **--t** _(1)_ - Number of threads to use.

## Output File
Each row is a gene and each column is a single string per sample per gene in the format “HAP_A_COUNT|HAP_B_COUNT”. Two output files are generated, one using only genome wide phased haplotype counts, such that the haplotype assignment is consistent across genes within an individual and with the input phased VCF. Another is generated that does not ensure genome wide haplotype phasing across genes, which includes more counts, but makes the haplotype assignment of A/B arbitrary and unrelated across genes within an individual or the VCF. Note: only the "gw_phased" output file may be used with phaser_cis_var.py as consistent phasing is required.

# phaser_cis_var.py
## Usage
Estimates effect sizes of regulatory variants using aggregated phASER haplotypic expression data. As input, takes a phASER haplotype expression matrix, a VCF, and a list of variants to calculate effect sizes for. To improve accuracy, the read-backed phased VCFs produced by phASER should be used, but first need to be combined across individuals, which can be performed using, e.g., “bcftools merge ind1.vcf.gz ind2.vcf.gz …”. The tool phases variant alleles with haplotype expression data and outputs numerous statistics, including allelic fold change (aFC, [Mohammadi, 2017](http://www.genome.org/cgi/doi/10.1101/gr.216747.116)) per sample, and a median across samples for individuals that heterozygous for the variant of interest. This median can be used as the estimate of the regulatory variant effect size. The output also includes statistics calculated for homozygous individuals, which can be used to test for increased allelic imbalance in heterozygotes as compared to homozygotes, as would be expected for true regulatory variants.

## Arguments - phaser_expr_matrix.py
### Required
* **--bed** - Outputted "gw_phased" expression matrix BED from phaser_expr_matrix.py
* **--vcf** - Phased VCF.
* **--pair** - File containing variant gene pairs to measure.
* **--map** - File containing VCF to BED sample mapping.
* **--o** - Output file.

### Optional
* **--pc** _(1)_ - Psuedocount to add to each haplotype count when calcualting aFC. Prevents division by zero errors.
* **--min_cov** _(8)_ - Minimum total coverage for a given individual's gene aFC to be included in calculations.
* **--chr** - Restrict to specific chromosome.
* **--bs** _(10000)_ - Number of bootstraps to generate 95% CIs.
* **--ignore_v** _(0)_ - Ignore version.
* **--t** _(1)_ - Number of threads to use.

## Output File

* 1 - **gene** - Gene ID
* 2 - **var_id** - Variant ID
* 3 - **var_chr** -  Variant chromosome
* 4 - **var_pos** - Variant position
* 5 - **var_het_n** - Number of heterozygous individuals
* 6 - **var_hom_n** - Number of homozygous individuals
* 7 - **het_hom_pvalue** - Ranksum test p-value comparing |aFCs| of heterozygous individuals vs homozygous individuals
* 8 - **var_het_afc_lower** - Lower 95% CI for heterozygous aFC
* 9 - **var_het_afc** - Median heterozygous aFC
* 10 - **var_het_afc_upper** - Upper 95% CI for heterozygous aFC
* 11 - **var_het_pval** - Two-sided empiral p-value against null that aFC overlaps zero
* 12 - **var_het_abs_afc_lower** - Lower 95% CI for heterozygous |aFC|
* 13 - **var_het_abs_afc** - Median heterozygous |aFC|
* 14 - **var_het_abs_afc_upper** - Upper 95% CI for heterozygous |aFC|
* 15 - **var_hom_afc_lower** - Lower 95% CI for homozygous aFC 
* 16 - **var_hom_afc** - Median homozygous aFC
* 17 - **var_hom_afc_upper** - Upper 95% CI for homozygous aFC 
* 18 - **var_hom_abs_afc_lower** - Lower 95% CI for homozygous |aFC| 
* 19 - **var_hom_abs_afc** - Median homozygous |aFC|
* 20 - **var_hom_abs_afc_upper** - Upper 95% CI for homozygous |aFC|
* 21 - **var_het_afcs** - Comma-seperated list of heterozygous individual calculated aFCs
* 22 - **var_hom_afcs** - Comma-seperated list of homozygous individual calculated aFCs
* 23 - **var_het_ref_counts** - Comma-seperated list of heterozygous individual reference variant haplotype counts
* 24 - **var_het_alt_counts** - Comma-seperated list of heterozygous individual alternative variant haplotype counts
* 25 - **var_hom_hap1_counts** - Comma-seperated list of homozygous individual haplotype 1 counts (note: distinction of haplotype 1/2 is arbitrary in variant homozygotes)
* 26 - **var_hom_hap2_counts** - Comma-seperated list of homozygous individual haplotype 2 counts (note: distinction of haplotype 1/2 is arbitrary in variant homozygotes)
* 27 - **var_het_sample_ids** - Comma-seperated list of individuals that are heterozygoys for the variant of interest
* 28 - **var_hom_sample_ids** - Comma-seperated list of individuals that are homozygous for the variant of interest

