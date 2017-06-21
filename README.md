# phASER
**ph**asing and **A**llele **S**pecific **E**xpression from **R**NA-seq

Performs haplotype phasing using read alignments in BAM format from both DNA and RNA based assays, and provides measures of haplotypic expression for RNA based assays.

Developed by [Stephane E. Castel](mailto:scastel@nygenome.org) in the [Lappalainen Lab](http://tllab.org) at the New York Genome Center and Columbia University Department of Systems Biology.

Please see our paper in [Nature Communications](http://www.nature.com/articles/ncomms12817) for details and benchmarking of the method.

phASER is made available under the [GNU GENERAL PUBLIC LICENSE v3](https://github.com/secastel/phaser/tree/master/LICENSE).

Documentation Pages: [phASER](https://github.com/secastel/phaser/tree/master/phaser), [phASER Annotate](https://github.com/secastel/phaser/tree/master/phaser_annotate), [phASER Gene AE](https://github.com/secastel/phaser/tree/master/phaser_gene_ae)

![alt tag](https://raw.github.com/secastel/phaser/master/docs/phaser_workflow.png)

# IMPORTANT NOTE - BUG FIX
A bug was introduced in version 0.9.8 (12/16/16) and fixed in version 0.9.9.4 (06/21/17) that caused problems with haplotypic counts when using the --haplo_count_blacklist argument. This bug affects the haplotypic counts generated (haplotypic_counts.txt), and any downstream analyses of those counts, including generating gene level haplotypic expression with phaser_gene_ae. If the --haplo_count_blacklist was not specified, then the results were not affected. In addition, a new "hg19_haplo_count_blacklist.bed.gz" file has been uploaded, which addresses problems related to this issue. If you used the --haplo_count_blacklist argument with a version of phASER between 0.9.8 and 0.9.9.3 you must re-run your analyses with version 0.9.9.4+.
