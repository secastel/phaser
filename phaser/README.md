
# phASER
**ph**asing and **A**llele **S**pecific **E**xpression from **R**NA-seq

Performs haplotype phasing using read alignments in BAM format from both DNA and RNA based assays.

Developed by [Stephane E. Castel](mailto:scastel@nygenome.org) in the [Lappalainen Lab](http://tllab.org) at the New York Genome Center and Columbia University Department of Systems Biology.

Runs on Python 2.7.x and has the following dependencies: [IntervalTree](https://github.com/jamescasbon/PyVCF), [pyVCF](https://pyvcf.readthedocs.org), [SciPy](http://www.scipy.org), [NumPy](http://www.numpy.org), [samtools](http://www.htslib.org), [bedtools](http://bedtools.readthedocs.org)

#Usage
Requires a BAM and VCF, produces a VCF with computed haplotype phases and a result file containing haplotype details, statistics, and read counts.

Test case: 1000 Genomes individual NA06986 with Phase 3 genotype data and Geuvadis RNA-Seq data.
* VCF: ftp://ftp.nygenome.org/sec/phaser/NA06986.vcf.gz
* BAM: ftp://ftp.nygenome.org/sec/phaser/NA06986.2.M_111215_4.bam

#Arguments
##Required
* **--bam** - Comma separated list of BAM files containing reads. Duplicates should be marked, and files should be indexed using samtools index.
* **--vcf** - VCF file containing genotype for the sample. May be gzipped. Chromosome names must be consistent between BAM and VCF.
* **--sample** - Name of sample to use in VCF file.
* **--baseq** - Minimum base quality at the SNP required for reads to be counted.
* **--mapq** - Mimimum mapping quality for reads to be counted. Can be a comma separated list, each value corresponding to the min MAPQ for a file in the input BAM list. Useful in cases when using both for example DNA and RNA libraries which having different mapping qualities.
* **--o** - Output file prefix name.

#Optional
* **--cc_threshold** _(0.01)_ - Threshold for significant conflicting variant configuration. The connection between any two variants with a conflicting configuration p-value lower than this threshold will be removed.
* **--isize** _(0)_ - Maximum allowed insert size for read pairs. Can be a comma separated list, each value corresponding to a max isize for a file in the input BAM list. Useful in cases when using both for example DNA and RNA libraries which having different expected insert sizes. Set to 0 for no maximum size.
* **--as_q_cutoff** _(0.05)_ - Bottom quantile to cutoff for alignment score. Reads with less than this alignment score quantile will not be included in the phasing.
* **--ab_q_cutoff** _(0)_ - Bottom quantile to cutoff for read aligned bases. Reads with fewer aligned bases than this aligned bases quantile will not be included in the phasing.
* **--blacklist** _()_ - BED file containing genomic intervals to be excluded from phasing (for example HLA).
* **--write_vcf** _(1)_ - Create a VCF containing phasing information (0,1).
* **--include_indels** _(0)_ - Include indels in the analysis (0,1). NOTE: since mapping is a problem for indels including them will likely result in poor quality phasing unless specific precautions have been taken.
* **--output_read_ids** _(0)_ - Output read IDs in the coverage files (0,1).
* **--output_orphans** _(0)_ - Output reads which fail to map to either allele (0, 1).
* **--remove_dups** _(1)_ - Remove duplicate reads from all analyses (0,1).
* **--pass_only** _(1)_ - Only use variants labled with PASS in the VCF filter field (0,1).
* **--min_cov** _(0)_ - Minimum total coverage level before outputting haplotypic counts.
* **--unphased_vars** _(1)_ - Output unphased variants (singletons) in the haplotypic_counts and haplotypes files (0,1). **NOTE** if you intend to run phASER Gene AE this must be enabled.

##Genome Wide Phasing
* **--gw_phase_method** _(0)_ - Method to use for determing genome wide phasing. NOTE requires input VCF to be phased, and optionally a VCF with allele frequencies (see --gw_af_vcf). 0 = Use most common haplotype phase. 1 = MAF weighted phase anchoring.
* **--gw_af_vcf** _()_ - VCF with allele frequencies from the population which was used to do the phasing. If left blank it will look for an allele frequency in the input VCF (--vcf).
* **--gw_af_field** _('AF')_ - Field from --gw_af_vcf to use for allele frequency.
* **--gw_phase_vcf** _(0)_ - Rephase GT field of output VCF using phASER genome wide phasing (0,1). See --gw_phase_method for options.
* **--gw_phase_vcf_min_confidence** _(0.90)_ - If replacing GT field in VCF only replace when phASER haplotype gw_confidence >= this value.

##Performance Related
* **--threads** _(1)_ - Maximum number of threads to use. Note the maximum thread count for some tasks is bounded by the data (for example 1 thread per contig for haplotype construction).
* **--max_block_size** _(15)_ - Maximum number of variants to phase at once. Number of haplotypes tested = 2 ^ # variants in block. Blocks larger than this will be split into sub blocks, phased, and then the best scoring sub blocks will be phased with each other.
* **--reads_mem** _(0)_ - Store reads that overlap variants in memory (0,1). NOTE: storing reads in memory increases speed, but grealtly increases memory overhead.
* **--vcf_mem** _(1)_ - Store the entire VCF into memory for quick processing (0,1). Should not be used with large, multisample VCF files.
* **--temp_dir** _()_ - Location of temporary directory to use for storing files. If left blank will default to system temp dir. NOTE: potentially large files will be stored in this directory, so please ensure there is sufficient free space.
* **--max_items_per_thread** _(100,000)_ - Maximum number of items that can be assigned to a single thread to process. NOTE: if this number is too high Python will stall when trying to join the pools.

##Debug / Development / Reporting
* **--show_warning** _(0)_ - Show warnings in stdout (0,1).
* **--debug** _(0)_ - Show debug mode messages (0,1).
* **--chr** _()_ - Restrict haplotype phasing to a specific chromosome.
* **--unique_ids** _(0)_ - Generate and output unique IDs instead of those provided in the VCF (0,1). NOTE: this should be used if your VCF does not contain a unique ID for each variant.
* **--output_network** _()_ - Output the haplotype connection network for the given variant.
* **--only_phased** _(0)_ - Only use variants in the input VCF which were marked as phased (0,1).

#Output Files

##*out_prefix*.haplotypes.txt

Contains all haplotypes phased along with detailed phasing statistics.

* 1 - **contig** - Contig haplotype is found on.
* 2 - **start** - Start position of haplotype (1 based).
* 3 - **stop** - Stop position of haplotype (1 based).
* 4 - **length** - Length of haplotype in bases (stop-start).
* 5 - **variants** - Number of variants phased in this haplotype.
* 6 - **variant_ids** -  Comma separated list of variant IDs phased in the haplotype. If --unique_ids enabled this will be the generated IDs.
* 7 - **variant_alleles** - Comma separated list of alleles found on each haplotype separated by "|".
* 8 - **reads_hap_a** - Number of unique reads mapping to haplotype A.
* 9 - **reads_hap_b** - Number of unique reads mapping to haplotype B.
* 10 - **reads_total** - Total number of unique reads covering this haplotype (reads_hap_a + reads_hap_b).
* 11 - **connections_supporting** - Number of variant connections (1 connection = 1 read spanning 2 variants) which support the chosen haplotype.
* 12 - **connections_total** - Total number of variant connections observed for this haplotype.
* 13 - **annotated_phase** - If input VCF contained phasing information, this will give the annotated phase for each haplotype separated by "|".
* 14 - **phase_concordant** - The phase is concordant with the input VCF annotation if each of the variants on a given haplotype the same annotated phase (in annotated_phase).
* 15 - **gw_phase** - Genome wide phasing for the alleles. The method for determing the genome wide phase is specified by --gw_phase_method. The phase will be outputted in the GT field of the VCF if --gw_phase_vcf = 1.
* 16 - **gw_confidence** - This represents the confidence of the GW phase assignment (between 0.5 and 1). A value of 0.5 indicates equal support for the two genome wide phasing configurations. In these cases the phase from the input VCF will be outputted.

##*out_prefix*.vcf.gz

Output of the input VCF file gzipped and containing the following added fields in the appropriate sample column:
* **GT** - This will only be outputted (overwriting the original input phase with PW) if --gw_phase_vcf = 1 and the phase statistics meet the criteria specified by --gw_phase_vcf_min_config_p and --gw_phase_vcf_min_confidence.
* **PW** - phASER Genome Wide Phase - The method for determing the genome wide phase is specified by --gw_phase_method. If no phASER phasing is available, or the there is no support for a specific phase (PC = 0.5) the input VCF phase will be outputted.
* **PC** - phASER Genome Wide Phase Confidence - This represents the confidence of the GW phase assignment (between 0.5 and 1). A value of 0.5 indicates equal support for the two genome wide phasing configurations. In these cases the phase from the input VCF will be outputted.
* **PG** - phASER Local Block Genotype - The phase of of this variant in the block specified in PB.
* **PI** - phASER Local Block Index - Unique index generated for each haplotype block.
* **PB** - phASER Local Block Variants - Comma separated list of variant IDs found in this haplotype block.

##*out_prefix*.allele_config.txt

This file contains the allele configuration for each variant where a phase could be established using the input reads. This file can be used to identify all cases of compound heterozygousity. Note: each variant pair will be listed twice, with A and B switched.

* **variant_a** - Unique variant ID for first variant. Format = chr_pos_ref_alt.
* **rsid_a** - RS ID for variant A.
* **variant_b** - Unique variant ID for second variant. Format = chr_pos_ref_alt.
* **rsid_b** - RS ID for variant B.
* **configuration** - Haplotype configuration for the two variants listed, one of two possibilities: "trans" = ref,alt|alt,ref (compound heterozygote), "cis" = ref,ref|alt,alt.

##*out_prefix*.allelic_counts.txt

Contains reference and alternative read counts for each heterozygous variant used for phasing. Format is the same as the [GATK ASEReadCounter](https://www.broadinstitute.org/gatk/gatkdocs/org_broadinstitute_gatk_tools_walkers_rnaseq_ASEReadCounter.php) and [allelecounter](https://github.com/secastel/allelecounter) outputs. If --min_cov is set only includes variants where totalCount ≥ min_cov.
* 1 - **contig** - Contig variant is found on.
* 2 - **position** - Position of variant, as per the VCF (1 based).
* 3 - **variantID** - ID of the variant. If --unique_ids enabled this will be the generated ID.
* 4 - **refAllele** - Reference allele base.
* 5 - **altAllele** - Alternate allele base.
* 6 - **refCount** - Reference allele read count.
* 7 - **altCount** - Alternate allele read count.
* 8 - **totalCount** - Total number of reads covering this base (refCount + altCount).

##*out_prefix*.haplotypic_counts.txt

Contains the number of unique reads that map to each haplotype for all phased haplotypes. If --min_cov is set only includes haplotypes where totalCount ≥ min_cov.

* 1 - **contig** - Contig haplotype is found on.
* 2 - **start** - Start position of haplotype (1 based).
* 3 - **stop** - Stop position of haplotype (1 based).
* 4 - **variants** - Comma separated list of variant IDs phased in the haplotype.
* 5 - **variantCount** - Number of variants contained in the haplotype.
* 6 - **haplotypeA** - Comma separated list of alleles found on haplotype A.
* 7 - **aCount** - Number of unique reads mapping to haplotype A.
* 8 - **haplotypeB** - Comma separated list of alleles found on haplotype B.
* 9 - **bCount** - Number of unique reads mapping to haplotype B.
* 10 - **totalCount** - Total number of unique reads covering this haplotype (aCount + bCount).

##*out_prefix*.variant_connections.txt

Statistics for every variant - variant connection observed by phASER in the data.

* **variant_a** - Unique ID of first variant.
* **variant_b** - Unique ID of second variant.
* **supporting_connections** - Number of reads which support the chosen phasing.
* **total_connections** - Total number of reads which overlap the two variants.
* **conflicting_configuration_p** - p value from the test for evidence for conflicting haplotype configuration.
* **phase_concordant** - If the input VCF was phased, is the phasing chosen by phASER concordant with the input phase (1/0).

##*out_prefix*.allele_config.txt

Configuration of all alleles with read backed phasing. Configuration can be either _cis_ or _trans_ and refers to the alternative alleles being either on the same haplotype (_cis_) or different haplotypes (_trans_).

* **variant_a** - Unique ID of first variant.
* **rsid_a** - RS ID of first variant.
* **variant_b** - Unique ID of second variant.
* **rsid_b** - RS ID of second variant.
* **configuration** - Haplotype configuration of alternative alleles.

##*out_prefix*.network.links.txt

If --output_network is enabled will contain the number of connections observed between each allele in the specific haplotype.

* **variantA** - VariantID:VariantAllele for first variant.
* **variantB** - VariantID:VariantAllele for second variant.
* **connections** - Number of connections observed between these two variants, 1 connection = 1 read spanning both variants.
* **inferred** - Number of inferred connections between these two variants, 1 inferred connection = 1 read spanning alternative alleles of two variants.

##*out_prefix*.network.nodes.txt

If --output_network is enabled will contain the names of each allele in the specific haplotype.

* **id** - VariantID:VariantAllele
