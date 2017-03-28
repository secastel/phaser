import pandas;
import gzip;
from intervaltree import IntervalTree;
import argparse;
import math;

def main():
	parser = argparse.ArgumentParser()
	# required
	parser.add_argument("--haplotypic_counts", required=True, help="Output file from phASER containing read counts for haplotype blocks. NOTE: unphased_vars must have been enabled when phASER was run.")
	parser.add_argument("--features", required=True, help="File in BED format (0 BASED COORDINATES - chr,start,stop,name) containing the features to produce counts for.")
	parser.add_argument("--o", required=True, help="Output file")
	
	# optional
	parser.add_argument("--gw_cutoff", type=float,default=0.9, help="Minimum genome wide phase confidence for phASER haplotype blocks.")
	parser.add_argument("--min_cov", type=int, default=0, help="Minimum total coverage for a feature to be outputted.")
	parser.add_argument("--no_gw_phase", type=int, default=0, help="Only use the haplotype block or SNP with maximum coverage per gene. Required if input VCF to phASER was unphased, --gw_cutoff will be ignored. NOTE with this option phasing between genes is not preserved (IE which haplotype is A/B is arbitrary and inconsistent between genes).")
	
	global args;
	args = parser.parse_args()
	
	version = "1.1.3";
	print("");
	print("##################################################")
	print("          Welcome to phASER Gene AE v%s"%(version));
	print("  Author: Stephane Castel (scastel@nygenome.org)")
	print("##################################################");
	print("");
	
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
		dict_feature_intervals[chrom][start:stop+1] = feature_index;
		
		dict_features[feature_index] = {'chr':chrom,'start':start,'stop':stop, 'name':name, 'aCount':0, 'bCount':0, 'variants':[]};
		
		feature_index += 1;
	
	print("#2 Loading haplotype counts...");
	# load haplotype counts
	# 	1 contig
	# 	2 start
	# 	3 stop
	# 	4 variants
	# 	5 variantCount
	# 	6 haplotypeA
	# 	7 haplotypeB
	# 	8 aCount
	# 	9 bCount
	# 	10 totalCount
	#	11 blockGWPhase
	#	12 gwStat

	df_haplo_counts = pandas.read_csv(args.haplotypic_counts, sep="\t", index_col=False);
	
	print("#3 Generating feature level haplotypic counts...");
	
	for index, row in df_haplo_counts.iterrows():
		# intersect haplotype with feature
		chrom = str(row['contig']);
		if row['totalCount'] > 0 and chrom in dict_feature_intervals:
			features = dict_feature_intervals[chrom][row['start']-1:row['stop']];
			
			if args.no_gw_phase == 1:
				for feature in features:
					feature_index = feature.data;
					if row['totalCount'] > (dict_features[feature_index]['aCount'] + dict_features[feature_index]['bCount']):
						dict_features[feature_index]['aCount'] = row['aCount'];
						dict_features[feature_index]['bCount'] = row['bCount'];
						dict_features[feature_index]['variants'] = row['variants'].split(",");
			
			elif row['blockGWPhase'] != "0/1" and float(row['gwStat'] >= args.gw_cutoff):
				for feature in features:
					feature_index = feature.data;
					if row['blockGWPhase'] == "0|1":
						# A = GW haplotype 0
						dict_features[feature_index]['aCount'] += row['aCount'];
						dict_features[feature_index]['bCount'] += row['bCount'];
					elif row['blockGWPhase'] == "1|0":
						# A = GW haplotype 1
						dict_features[feature_index]['aCount'] += row['bCount'];
						dict_features[feature_index]['bCount'] += row['aCount'];
				
					dict_features[feature_index]['variants'] += row['variants'].split(",");
	
	print("#4 Outputting feature haplotype counts...");	
	
	stream_out = open(args.o, "w");
	stream_out.write("\t".join(["contig","start","stop","name","aCount","bCount","totalCount","log2_aFC","n_variants","variants"])+"\n");
	
	for index in range(0,len(dict_features.keys())):
		total_cov = dict_features[index]['aCount']+dict_features[index]['bCount'];
		n_variants = len(dict_features[index]['variants']);
		log2_afc = zero_log(zero_divide(dict_features[index]['aCount'],dict_features[index]['bCount']),2);
		
		if total_cov >= args.min_cov:
			stream_out.write("\t".join(map(str,[dict_features[index]['chr'],dict_features[index]['start'],dict_features[index]['stop'],dict_features[index]['name'],dict_features[index]['aCount'],dict_features[index]['bCount'],total_cov, log2_afc, n_variants, ",".join(dict_features[index]['variants'])]))+"\n");
		else:
			stream_out.write("\t".join(map(str,[dict_features[index]['chr'],dict_features[index]['start'],dict_features[index]['stop'],dict_features[index]['name'],0,0,0, float('nan'), 0, ""]))+"\n");
			
	stream_out.close();

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
	
