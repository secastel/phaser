# phASER
Performs haplotype phasing using read alignments in BAM format from both DNA and RNA based experiments.

Developed by [Stephane E. Castel](mailto:scastel@nygenome.org) at the New York Genome Center and Columbia University Department of Systems Biology in the [Lappalainen Lab](http://tlab.org).

Runs on Python 2.7.x and has the following dependencies: [IntervalTree](https://github.com/jamescasbon/PyVCF), [SciPy](http://www.scipy.org), [NumPy](http://www.numpy.org), [Samtools](http://www.htslib.org), [Bedtools](http://bedtools.readthedocs.org)

#Usage
Requires a BAM and VCF, produces a VCF with computed haplotype phases and a result file containing haplotype details, statistics, and read counts.

#Arguments
##Required
* **--bam** - BAM file containing reads. Duplicates should be marked, and file should be indexed using Samtools index.
* **--vcf** - VCF file containing genotype for the sample. Must be gzipped and Tabix indexed. Chromosome names must be the same between BAM and VCF.
* **--sample** - Name of sample to use in VCF file.
* **--baseq** - Minimum base quality at the SNP required for reads to be counted.
* **--mapq** - Mimimum mapping qualityfor reads to be counted.
* **--o** - Output file prefix name.

#Optional
* **--blacklist** _()_ - BED file containing genomic intervals to be excluded from phasing (for example HLA).
* **--write_vcf** _(1)_ - Create a VCF containing phasing information (0,1).
* **--include_indels** _(0)_ - Include indels in the analysis (0,1). NOTE: since mapping is a problem for indels including them will likely result in poor quality phasing unless specific precautions have been taken.
* **--output_read_ids** _(0)_ - Output read IDs in the coverage files (0,1).
* **--output_orphans** _(0)_ - Output reads which fail to map to either allele (0, 1).
* **--remove_dups** _(1)_ - Remove duplicate reads from all analyses (0,1).
* **--pass_only** _(1)_ - Only use variants labled with PASS in the VCF filter field (0,1).
* **--min_cov** _(0)_ - Minimum total coverage level before outputting haplotypic counts.

##Performance Related
* **--threads** _(1)_ - Maximum number of threads to use. Note the maximum thread count for some tasks is bounded by the data (for example 1 thread per contig for haplotype construction).
* **--max_block_size** _(15)_ - Maximum number of variants to phase at once. Number of haplotypes tested = 2 ^ # variants in block. Blocks larger than this will be split into sub blocks, phased, and then the best scoring sub blocks will be phased with each other.
* **--reads_mem** _(0)_ - Store reads that overlap variants in memory (0,1). NOTE: storing reads in memory increases speed, but grealtly increases memory overhead.
* **--vcf_mem** _(1)_ - Store the entire VCF into memory for quick processing (0,1). Should not be used with large, multisampling VCF files.
* **--temp_dir** _()_ - Location of temporary directory to use for storing files. If left blank will default to system temp dir. NOTE: potentially large files will be stored in this directory, so please ensure there is sufficient free space.
* **--max_items_per_thread** _(100,000)_ - Maximum number of items that can be assigned to a single thread to process. NOTE: if this number is too high Python will stall when trying to join the pools.

##Debug / development / reporting
* **--show_warning** _(0)_ - Show warnings in stdout (0,1).
* **--debug** _(0)_ - Show debug mode messages (0,1).
* **--chr** _()_ - Restrict haplotype phasing to a specific chromosome.
* **--unique_ids** _(0)_ - Generate and output unique IDs instead of those provided in the VCF (0,1). NOTE: this should be used if your VCF does not contain a unique ID for each variant.
* **--output_network** _()_ - Output the haplotype connection network for the given variant.
* **--only_phased** _(0)_ - Only use variants in the input VCF which were marked as phased (0,1).

#Output Files
