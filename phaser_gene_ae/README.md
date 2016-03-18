# phASER Gene AE
Uses output from phASER to produce gene level haplotype counts for allelic expression studies. It does this by summing reads from both single variants and phASER haplotype blocks using their phase for each gene.

Developed by [Stephane E. Castel](mailto:scastel@nygenome.org) in the [Lappalainen Lab](http://tllab.org) at the New York Genome Center and Columbia University Department of Systems Biology.

Runs on Python 2.7.x and has the following dependencies: [pandas](http://pandas.pydata.org), [IntervalTree](https://github.com/jamescasbon/PyVCF)

#Usage
Requires phASER to have been run with a phased VCF as input with unphased_vars enabled. Takes an input BED format file containing the coordinates for genes (feautres) where haplotypic counts are to be measured.

**Important Note** - By default, the haplotypic counts produced by phASER are summed across all of the input libraries. This means that if you used, for example, both DNA and RNA input libraries, the counts produced in o.haplotypic_counts.txt would not be useful for allelic expression studies. In such cases the "--haplo_count_bam" argument can be used to specify the libraries that should be used to generate haplotypic counts.

#Arguments
##Required
* **--haplotypic_counts** - Output file from phASER containing read counts for haplotype blocks.
* **--features** - File in BED format (0 BASED COORDINATES - chr,start,stop,name) containing the features to produce counts for.
* **--o** - Output file.

##Optional
* **--min_cov** _(1)_ - Minimum total coverage for a feature to be outputted.
* **--gw_cutoff** _(0.9)_ - Minimum genome wide phase confidence for phASER haplotype blocks.
* **--no_gw_phase** _(0)_ - Only use the haplotype block or SNP with maximum coverage per gene. Required if input VCF to phASER was unphased, --gw_cutoff will be ignored. NOTE with this option phasing between genes is not preserved (IE which haplotype A/B is arbitrary and inconsistent between genes).

#Output File

Contains the haplotype counts (A = genome wide haplotype 0, B = genome wide haplotype 1) for each feature. Note that because global haplotypes are used counts between features can even be compared. For example for two features on the same chromosome, reads from haplotype A for each would have come from the same DNA molecule.

* 1 - **contig** - Feature contig.
* 2 - **start** - Feature start (0 based).
* 3 - **stop** - Feature stop (0 based).
* 4 - **name** - Feature name.
* 5 - **aCount** - Total allelic count for haplotype A.
* 6 - **bCount** - Total allelic count for haplotype B.
* 7 - **totalCount** - Total allelic coverage of this feature (aCount + bCount).
* 8 - **n_variants** - Number of variants with allelic data in this feature.
* 9 - **variants** - List of variants with allelic data in this feature (contig_position_ref_alt).
