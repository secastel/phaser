# The code below shows how to use phaser_cis_var.py to calculate aFCs from phASER haplotypic data for GTEx v8 Whole Blood top eQTLs

# Required File - Download Location
# phASER_WASP_GTEx_v8_matrix.gw_phased.txt.gz - GTEx Portal (https://gtexportal.org/home/datasets):
# phASER_GTEx_v8_merged.vcf.gz - dbGaP (phs000424)

python phaser_cis_var.py --bed phASER_WASP_GTEx_v8_matrix.gw_phased.txt.gz --vcf phASER_GTEx_v8_merged.vcf.gz --pairs Whole_Blood.test_pairs.txt --map Whole_Blood.sample_map.txt --o Whole_Blood.results.txt --ignore_v 1