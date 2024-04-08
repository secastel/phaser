import pandas;
import gzip;
from intervaltree import IntervalTree;
import argparse;
import math;
import sys;

def main():
	parser = argparse.ArgumentParser()
	# required
	parser.add_argument("--haplotypic_counts", required=True, help="Output file from phASER containing read counts for haplotype blocks. NOTE: unphased_vars must have been enabled when phASER was run.")
	parser.add_argument("--features", required=True, help="File in BED format (0 BASED COORDINATES - chr,start,stop,name) containing the features to produce counts for.")
	parser.add_argument("--o", required=True, help="Output file")

	# optional
	parser.add_argument("--id_separator", default="_", help="Separator used for generating unique variant IDs when phASER was run.")
	parser.add_argument("--gw_cutoff", type=float,default=0.9, help="Minimum genome wide phase confidence for phASER haplotype blocks.")
	parser.add_argument("--min_cov", type=int, default=0, help="Minimum total coverage for a feature to be outputted.")
	parser.add_argument("--min_haplo_maf", type=float, default=0, help="The minimum MAF used to phase a haplotype for it to be considered genome wide phased when generating gene level counts. Setting this number higher will result in more confident phasing if genotypes were population prephased. Value must be between 0 and 0.5.")

	global args;
	args = parser.parse_args()

	version = "1.2.0";
	print("");
	print("##################################################")
	print(("          Welcome to phASER Gene AE v%s"%(version)));
	print("  Author: Stephane Castel (scastel@nygenome.org)")
	print("##################################################");
	print("");

	if args.min_haplo_maf < 0 or args.min_haplo_maf > 0.5:
		print("ERROR - invalid value for min_haplo_maf specified. Value must be between 0 and 0.5.")
		sys.exit(1)

	global dict_features;
	print("#1 Loading features...");
	dict_feature_intervals = {};
	dict_features = {};
	stream_features = open(args.features, "r")

	feature_index = 0;
	for line in stream_features:
		columns = line.rstrip().split("\t");
		chrom = columns[0];
		start = int(columns[1]);
		stop = int(columns[2]);
		name = columns[3];

		if chrom not in dict_feature_intervals: dict_feature_intervals[chrom] = IntervalTree();
		dict_feature_intervals[chrom][start:stop] = feature_index;

		dict_features[feature_index] = {'chr':chrom,'start':start,'stop':stop, 'name':name, 'aCount':0, 'bCount':0, 'unphased_aCount':0,'unphased_bCount':0,'unphased_variants':[],'variants':[]};

		feature_index += 1;

	print("#2 Loading haplotype counts...");
	# 1 contig
	# 2 start
	# 3 stop
	# 4 variants
	# 5 variantCount
	# 6 variantsBlacklisted
	# 7 variantCountBlacklisted
	# 8 haplotypeA
	# 9 haplotypeB
	# 10 aCount
	# 11 bCount
	# 12 totalCount
	# 13 blockGWPhase
	# 14 gwStat
	# 15 max_haplo_maf
	# 16 bam
	# 17 aReads
	# 18 bReads

	df_haplo_counts_master = pandas.read_csv(args.haplotypic_counts, sep="\t", index_col=False);

	if "bam" not in df_haplo_counts_master.columns:
		print("ERROR - this version of phaser_gene_ae is only compatible with results from phASER v1.0.0+");
		sys.exit(1)

	stream_out = open(args.o, "w");
	stream_out.write("\t".join(["contig","start","stop","name","aCount","bCount","totalCount","log2_aFC","n_variants","variants","gw_phased","bam"])+"\n");

	print("#3 Processing results...")
	# produce a separate output file for each bam
	for xbam in set(df_haplo_counts_master['bam']):
		print(("    BAM: %s"%(xbam)))
		print("          generating feature level haplotypic counts...");
		maf_filtered = 0;
		df_haplo_counts = df_haplo_counts_master[(df_haplo_counts_master.bam == xbam)];
		# reset counts on features to 0
		for xfeature in dict_features:
			dict_features[xfeature]['aCount'] = 0;
			dict_features[xfeature]['bCount'] = 0;
			dict_features[xfeature]['variants'] = [];
			dict_features[xfeature]['unphased_aCount'] = 0;
			dict_features[xfeature]['unphased_bCount'] = 0;
			dict_features[xfeature]['unphased_variants'] = "";

		for index, row in df_haplo_counts.iterrows():
			# intersect haplotype with feature
			chrom = str(row['contig']);
			if row['totalCount'] > 0 and chrom in dict_feature_intervals:
				# BED coordinates are 0 based, haplotypes are 1 based
				features = dict_feature_intervals[chrom][row['start']-1:row['stop']];

				for feature in features:
					feature_index = feature.data;
					mapped_reads = variant_feature_reads(row,feature);

					if row['blockGWPhase'] != "0/1" and float(row['gwStat'] >= args.gw_cutoff):
						# if haplotype mafs have been provided only use haplotypes where the maximum MAF is greater than args.min_haplo_maf
						if args.min_haplo_maf > 0 and "max_haplo_maf" in df_haplo_counts.columns and row['max_haplo_maf'] < args.min_haplo_maf:
							# treat as unphased
							if mapped_reads['totalCount'] > (dict_features[feature_index]['unphased_aCount'] + dict_features[feature_index]['unphased_bCount']):
								dict_features[feature_index]['unphased_aCount'] = mapped_reads['aCount'];
								dict_features[feature_index]['unphased_bCount'] = mapped_reads['bCount'];
								dict_features[feature_index]['unphased_variants'] = mapped_reads['variants'];
							maf_filtered += 1;
							continue;

						if row['blockGWPhase'] == "0|1":
							# A = GW haplotype 0
							dict_features[feature_index]['aCount'] += mapped_reads['aCount'];
							dict_features[feature_index]['bCount'] += mapped_reads['bCount'];
						elif row['blockGWPhase'] == "1|0":
							# A = GW haplotype 1
							dict_features[feature_index]['aCount'] += mapped_reads['bCount'];
							dict_features[feature_index]['bCount'] += mapped_reads['aCount'];

						dict_features[feature_index]['variants'] += mapped_reads['variants'];

					else:
						# this is an unphased block
						if mapped_reads['totalCount'] > (dict_features[feature_index]['unphased_aCount'] + dict_features[feature_index]['unphased_bCount']):
							dict_features[feature_index]['unphased_aCount'] = mapped_reads['aCount'];
							dict_features[feature_index]['unphased_bCount'] = mapped_reads['bCount'];
							dict_features[feature_index]['unphased_variants'] = mapped_reads['variants'];

		if maf_filtered > 0:
			print(("          %d of %d haplotypes treated as unphased due to low MAF"%(maf_filtered, len(df_haplo_counts.index))));

		print("          outputting feature haplotype counts...");

		for index in range(0,len(list(dict_features.keys()))):
			# decide whether to use unphased or phased counts
			# use whichever is higher
			if (dict_features[index]['aCount'] + dict_features[index]['bCount']) >= (dict_features[index]['unphased_aCount'] + dict_features[index]['unphased_bCount']):
				# use phased counts
				total_cov = dict_features[index]['aCount']+dict_features[index]['bCount'];
				n_variants = len(dict_features[index]['variants']);
				log2_afc = zero_log(zero_divide(dict_features[index]['aCount'],dict_features[index]['bCount']),2);

				if total_cov >= args.min_cov:
					stream_out.write("\t".join(map(str,[dict_features[index]['chr'],dict_features[index]['start'],dict_features[index]['stop'],dict_features[index]['name'],dict_features[index]['aCount'],dict_features[index]['bCount'],total_cov, log2_afc, n_variants, ",".join(dict_features[index]['variants']),1,xbam]))+"\n");
			elif (dict_features[index]['aCount'] + dict_features[index]['bCount']) < (dict_features[index]['unphased_aCount'] + dict_features[index]['unphased_bCount']):
				# use unphased count (this is just the highest covered haplotype block or single variant)
				total_cov = dict_features[index]['unphased_aCount']+dict_features[index]['unphased_bCount'];
				n_variants = len(dict_features[index]['unphased_variants']);
				log2_afc = zero_log(zero_divide(dict_features[index]['unphased_aCount'],dict_features[index]['unphased_bCount']),2);
				if total_cov >= args.min_cov:
					stream_out.write("\t".join(map(str,[dict_features[index]['chr'],dict_features[index]['start'],dict_features[index]['stop'],dict_features[index]['name'],dict_features[index]['unphased_aCount'],dict_features[index]['unphased_bCount'],total_cov, log2_afc, n_variants, ",".join(dict_features[index]['unphased_variants']),0,xbam]))+"\n");
			elif args.min_cov == 0:
				stream_out.write("\t".join(map(str,[dict_features[index]['chr'],dict_features[index]['start'],dict_features[index]['stop'],dict_features[index]['name'],0,0,0, float('nan'), 0, "",float('nan'),xbam]))+"\n");


	stream_out.close();


def variant_feature_reads(row,feature):
	global dict_features;
	# unpack haplotype reads, only count reads from variants that overlap the feature
	hap_a_reads = [];
	hap_b_reads = [];
	used_vars = [];

	xvars = row['variants'].split(",");

	# check that the correct ID separator has been provided
	if args.id_separator not in xvars[0] or xvars[0].count(args.id_separator) < 3:
		print("ERROR - ID separator not found in variant ID, please ensure that --id_separator is set correctly.")
		sys.exit(1)

	# find those variants that overlap the feature, record their reads
	for xvar in xvars:
		xvar_index = xvars.index(xvar);

		fields = xvar.split(args.id_separator);
		xvar_pos = int(fields[1]);

		# check if the variant overlaps the feature
		if (xvar_pos-1) - feature.begin >= 0 and (xvar_pos-1) - feature.end <= 0:
			# if so add the reads from this variant
			used_vars.append(xvar);

			if len(xvars) == 1:
				# if haplotype only has one variant then there are no read names since they are not needed
				# generate some fake read ids
				hap_a_reads += list(map(str,list(range(0,int(row['aCount'])))));
				hap_b_reads += list(map(str,list(range(0,int(row['bCount'])))));
			else:
				hap_a_reads += str(row['aReads']).split(";")[xvar_index].split(",");
				hap_b_reads += str(row['bReads']).split(";")[xvar_index].split(",");

	# make set of read ids, results in no double counting of reads / molecules
	hap_a_reads = set(hap_a_reads);
	hap_b_reads = set(hap_b_reads);

	# remove blank read IDs (these are created when there are no reads mapping to a variant on a given haplotype)
	if "" in hap_a_reads: hap_a_reads.remove("");
	if "" in hap_b_reads: hap_b_reads.remove("");

	# use a set so that reads are not double counted if they overlap more than one variant
	aCount = len(set(hap_a_reads));
	bCount = len(set(hap_b_reads));

	return({'variants':used_vars,'aCount':aCount,'bCount':bCount, 'totalCount':aCount+bCount});

def safe_dict(key,xdict,null_val):
	if key in xdict:
		return(xdict[key]);
	else:
		return(null_val);

def zero_divide(a,b):
	if b == 0:
		return(float('inf'));
	else:
		return(float(a)/float(b));

def zero_log(value,base):
	if value == 0:
		return(float('-inf'));
	else:
		return(math.log(value,base));

if __name__ == "__main__":
	main();
