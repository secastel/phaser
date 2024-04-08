# phASER Gene AE
Uses output from phASER to produce gene level haplotype counts for allelic expression studies. It does this by summing reads from both single variants and phASER haplotype blocks using their phase for each gene.

Developed by [Stephane E. Castel](mailto:stephanecastel@gmail.com) in the [Lappalainen Lab](http://tllab.org) at the New York Genome Center and Columbia University Department of Systems Biology.

Runs on Python 3.x.x and has the following dependencies: [pandas](http://pandas.pydata.org), [IntervalTree](https://github.com/jamescasbon/PyVCF)

# Usage
Requires phASER to have been run with a phased VCF as input with unphased_vars enabled. Takes an input BED format file containing the coordinates for genes (feautres) where haplotypic counts are to be measured. **NOTE** this version of phaser_gene_ae is only compatible with results from phASER v1.0.0+.

**Important Note** - If multiple input BAMs were used when running phASER, gene level haplotypic counts will be generated **for each input BAM independently**. In the output file, the column "bam" indicates which input BAM the haplotypic counts correspond to.

**Useful files**

The specific features to produce haplotypic counts for must be provided in BED format. This is most often genes. A file containing coordinates for gencode genes is included here for convenience. Note that these annotations have been filtered to be consistent with what was used for calling GTEx eQTLs. Please ensure that they contain annotations appropriate for your specific analysis.

hg19:
* Without 'chr' in contig name: https://www.dropbox.com/s/1u9zo1kx61zx6ca/gencode.v19.GRCh37.genes.bed.gz?dl=0
* With 'chr' in contig name: https://www.dropbox.com/s/am09zwpjhs01k8u/gencode.v19.GRCh37.genes.chr.bed.gz?dl=0

hg38:
* Without 'chr' in contig name: https://www.dropbox.com/s/ncdika74bcgtrc2/gencode.v26.GRCh38.genes.bed.gz?dl=0
* With 'chr' in contig name: https://www.dropbox.com/s/4vmgbk9fcegxk0r/gencode.v26.GRCh38.genes.chr.bed.gz?dl=0


# Arguments
## Required
* **--haplotypic_counts** - Output file from phASER containing read counts for haplotype blocks.
* **--features** - File in BED format (0 BASED COORDINATES - chr,start,stop,name) containing the features to produce counts for.
* **--o** - Output file.

## Optional
* **--id_separator** _(\_)_ - Separator used for generating unique variant IDs when phASER was run.
* **--gw_cutoff** _(0.9)_ - Minimum genome wide phase confidence for phASER haplotype blocks.
* **--min_cov** _(0)_ - Minimum total coverage for a feature to be outputted.
* **--min_haplo_maf** _(0)_ - The minimum MAF used to phase a haplotype for it to be considered genome wide phased when generating gene level counts. Setting this number higher will result in more confident phasing if genotypes were population prephased. Value must be between 0 and 0.5.

# Output File

Contains the haplotype counts (A = genome wide haplotype 0, B = genome wide haplotype 1) for each feature, for each input BAM. The column "gw_phased" indicates if a feature is genome wide phased, meaning that for two features that are both genome wide phased, reads from haplotype A for each would have come from the same DNA molecule.

* 1 - **contig** - Feature contig.
* 2 - **start** - Feature start (0 based).
* 3 - **stop** - Feature stop (0 based).
* 4 - **name** - Feature name.
* 5 - **aCount** - Total allelic count for haplotype A.
* 6 - **bCount** - Total allelic count for haplotype B.
* 7 - **totalCount** - Total allelic coverage of this feature (aCount + bCount).
* 8 - **log2_aFC** - Effect size for the allelic imbalance reported as allelic fold change (log2(aCount/bCount)) defined in our [paper](http://biorxiv.org/content/early/2016/09/30/078717).
* 9 - **n_variants** - Number of variants with allelic data in this feature.
* 10 - **variants** - List of variants with allelic data in this feature (contig_position_ref_alt).
* 11 - **gw_phased** - The feature is genome wide phased, meaning that for two features that are both genome wide phased, reads from haplotype A for each would have come from the same DNA molecule (0,1).
* 12 - **bam** - Input BAM used when running phASER that these counts correspond to.
