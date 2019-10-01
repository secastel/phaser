import pandas;
import argparse;
import subprocess;
from os.path import isfile;
import numpy;
import multiprocessing;
import pysam;
from scipy.stats import ranksums;
from math import log;
import tempfile;
from pandas.compat import StringIO
import gzip;

def main():
	parser = argparse.ArgumentParser();
	## required
	parser.add_argument("--bed", type=str, required=True, help="phASER gw_phased expression matrix BED");
	parser.add_argument("--vcf", type=str, required=True, help="Phased VCF");
	parser.add_argument("--pairs", type=str, required=True, help="File containing variant gene pairs to measure");
	parser.add_argument("--map", type=str, required=True, help="File containing VCF to BED sample mapping");
	parser.add_argument("--o", type=str, required=True, help="Output");
	
	## optional
	parser.add_argument("--pc", default=1, type=int, help="Psuedocount to add to each haplotype count when calcualting aFC. Prevents division by zero errors.")
	parser.add_argument("--min_cov", type=int, default=8, help="Minimum total coverage for a given individual's gene aFC to be included in calculations.")
	parser.add_argument("--chr", type=str, required=False, default="", help="Restrict to specific chromosome");
	parser.add_argument("--bs", type=int, required=False, default=10000, help="Number of bootstraps to generate 95% CIs");
	parser.add_argument("--ignore_v", type=int, required=False, default=0, help="Ignore gene version");
	parser.add_argument("--t", type=int, required=False, default=1, help="Threads");
	
	version = "0.1.0";
	print("");
	print("##################################################")
	print("          Welcome to phASER-POP v%s"%(version));
	print("  Author: Stephane Castel (scastel@nygenome.org)")
	print("##################################################");
	print("");

	global args;
	args = parser.parse_args();

	print("#1 Loading sample map...")
	df_map = pandas.read_csv(args.map, sep="\t", index_col=False);
	global dict_map;
	dict_map = dict(zip(df_map['vcf_sample'],df_map['bed_sample']));

	print("#2 Loading variant gene pairs...")
	global df_pairs;
	df_pairs = pandas.read_csv(args.pairs, sep="\t", index_col=False);
	if args.ignore_v == 1: df_pairs['gene_id'] = [x.split(".")[0] for x in df_pairs['gene_id']];
	if args.chr != "":
		df_pairs['var_contig'] = map(str, df_pairs['var_contig'])
		df_pairs = df_pairs[(df_pairs.var_contig == args.chr)];

	print("#3 Loading phASER BED...")
	global df_phaser;
	in_bed = args.bed;
	if args.chr != "":
		print("     subsetting chr %s from input BED..."%(args.chr));
		xfile = tempfile.NamedTemporaryFile(delete=False);
		in_bed = xfile.name;
		xfile.close();
		subprocess.call("tabix -h "+args.bed+" "+args.chr+": > "+in_bed, shell=True)

	# filter the expression matrix to only load those genes we actually need to do the analysis
	xfile = tempfile.NamedTemporaryFile(delete=False);
	if ".gz" in in_bed:
		stream_in = gzip.open(in_bed, "r");
	else:
		stream_in = open(in_bed, "r");

	use_lines = [];
	set_use_genes = set(df_pairs['gene_id'].tolist());
	for xline in stream_in:
		xline = xline.rstrip();
		if xline.startswith("#"):
			use_lines.append(xline);
		xcols = xline.split("\t");
		if args.ignore_v == True: xcols[1] = xcols[1].split(".")[0];
		if xcols[1] in set_use_genes:
			use_lines.append(xline);

	df_phaser = pandas.read_csv(StringIO("\n".join(use_lines)), sep="\t");
	df_phaser.index = df_phaser['name'];
	if args.ignore_v == 1: df_phaser.index = [x.split(".")[0] for x in df_phaser.index];

	if len(df_phaser.index) == 0:
		print("     ERROR: no phASER data read from input... check that chromsome naming is correct.")
		quit();

	print("#4 Measuring allelic expression...")
	xpool = chunks(df_pairs.index, len(df_pairs.index)/args.t);
	xresults = parallelize(measure_effect, [x for x in xpool]);

	list_results = [];
	for xresult in xresults:
		for xitem in xresult:
			list_results.append(xitem);

	df_result = pandas.DataFrame(list_results, columns=['xindex','gene','var_id','var_chr','var_pos','var_het_n','var_hom_n','het_hom_pvalue','var_het_afc_lower','var_het_afc','var_het_afc_upper','var_het_pval','var_het_abs_afc_lower','var_het_abs_afc','var_het_abs_afc_upper','var_hom_afc_lower','var_hom_afc','var_hom_afc_upper','var_hom_abs_afc_lower','var_hom_abs_afc','var_hom_abs_afc_upper','var_het_afcs','var_hom_afcs','var_het_ref_counts','var_het_alt_counts','var_hom_hap1_counts','var_hom_hap2_counts','var_het_sample_ids','var_hom_sample_ids']);
	
	# sort so that output order is the same as input
	df_result = df_result.sort_values(by=['xindex']);
	df_result = df_result.drop(['xindex'], axis=1);

	df_result.to_csv(args.o, sep="\t", index=False);

	if args.chr != "":
		subprocess.call("rm "+in_bed, shell=True)

def measure_effect(xindices):
	global args;
	global dict_map;
	global df_phaser;
	global df_pairs;

	out_result = [];

	# initialize VCF reader for retrieving genotypes
	xvcf_reader = vcf_reader(args.vcf);

	for xindex in xindices:
		row_test = df_pairs.loc[xindex];
		if row_test['gene_id'] in df_phaser.index:
			afcs = [[],[]];
			phaser_counts = [[[],[]],[[],[]]];
			ids = [[],[]];
			row_phaser = df_phaser.loc[row_test['gene_id']];
			# retrieve genotypes from VCF

			records = xvcf_reader.retrieve_variant(row_test['var_contig'],row_test['var_pos']);
			for xrecord in records:
				# if var_ref and var_alt have been set match on those otherwise try on ID
				if (row_test['var_ref'] != "" and row_test['var_alt'] != ""  and xrecord['REF'] == row_test['var_ref'] and xrecord['ALT'] == row_test['var_alt']) or (xrecord['ID'] == row_test['var_id']):
					gt_index = xrecord['FORMAT'].split(":").index("GT");
					# got the variant, now retrieve each individual's genotype
					for xsamp in dict_map.keys():
						if xsamp in xrecord and dict_map[xsamp] in df_phaser.columns:
							xgt = xrecord[xsamp].split(":")[gt_index];
							# need phased genotype data
							if "|" in xgt:
								counts = map(float,row_phaser[dict_map[xsamp]].split("|"));
								if sum(counts) >= args.min_cov:
									afc_phaser = log(float(counts[0]+args.pc)/float(counts[1]+args.pc),2);
									if "0" in xgt and "1" in xgt:
										# phASER aFCs are listed as hap0/hap1
										# we want aFC relative to the SNP of interest such that it is ALT/REF
										# so if the SNP is on hap1 then flip the aFC
										alt_index = xgt.split("|").index("1");
										if alt_index == 1:afc_phaser *= -1;
										afcs[0].append(afc_phaser);
										ids[0].append(xsamp)
										phaser_counts[0][0].append(int(counts[int(not alt_index)]));
										phaser_counts[0][1].append(int(counts[alt_index]));
									elif xgt.count("0") == 2 or xgt.count("1") == 2:
										afcs[1].append(afc_phaser);
										ids[1].append(xsamp)
										phaser_counts[1][0].append(int(counts[0]));
										phaser_counts[1][1].append(int(counts[1]));
					
					# calculate |aFC|
					abs_afcs = [[],[]];
					abs_afcs[0] = map(abs,afcs[0]);
					abs_afcs[1] = map(abs,afcs[1]);

					# bootstrap 95% CIs
					het_ci = bootstrap_ci(afcs[0],args.bs,numpy.median,True)
					het_abs_ci = bootstrap_ci(abs_afcs[0],args.bs,numpy.median)

					hom_ci = bootstrap_ci(afcs[1],args.bs,numpy.median)
					hom_abs_ci = bootstrap_ci(abs_afcs[1],args.bs,numpy.median)			
					
					stat,pval = ranksums(abs_afcs[0],abs_afcs[1]);

					out_result.append([xindex,row_phaser['name'],row_test['var_id'],row_test['var_contig'],row_test['var_pos'],len(afcs[0]),len(afcs[1]),pval]+het_ci+het_abs_ci+hom_ci+hom_abs_ci+[list_to_str(afcs[0]),list_to_str(afcs[1]),list_to_str(phaser_counts[0][0]),list_to_str(phaser_counts[0][1]),list_to_str(phaser_counts[1][0]),list_to_str(phaser_counts[1][1]),list_to_str(ids[0]),list_to_str(ids[1])]);

	return(out_result);

def chunks(l, n):
	for i in range(0, len(l), n):
		yield l[i:i + n];

def parallelize(function, pool_input):
	global args;

	if len(pool_input) > 0:
		threads = min([len(pool_input),args.t]);
		if args.t > 1:
			pool = multiprocessing.Pool(processes=threads);
			pool_output = pool.map(function, pool_input);
			pool.close() # no more tasks
			pool.join()  # wrap up current tasks
		else:
			pool_output = [];
			for input in pool_input:
				pool_output.append(function(input));
	else:
		pool_output = [];

	return(pool_output);

def bootstrap_ci(xinput,xbs,xfun,return_p = False):
	if len(xinput) > 0:
		vals = [];
		for i in range(0,xbs):
			vals.append(xfun(numpy.random.choice(xinput,replace=True,size=len(xinput))))

		# calculate_pvalue
		xout = [numpy.percentile(vals,2.5),xfun(xinput),numpy.percentile(vals, 97.5)];
		if return_p == True:
			p_val = (float(min([sum([int(x>0) for x in vals]),sum([int(x<0) for x in vals])]))/float(xbs)) * 2;
			xout = xout + [p_val];

		return(xout);
	else:
		xout = [float('nan'),float('nan'),float('nan')];
		if return_p == True: xout = xout + [float('nan')];
		
		return(xout);

def list_to_str(xlist, sep=","):
	return(sep.join(map(str,xlist)));

class vcf_reader:
	def __init__(self, vcf_path):
		# setup a reader
		self.tabix_vcf = pysam.Tabixfile(vcf_path,"r");
		# get the column names
		col_names = subprocess.check_output("tabix -H "+vcf_path+" | grep 'CHROM'", shell=True).replace("#","").rstrip().split("\t");
		self.dict_cols = dict(zip(range(0,len(col_names)),col_names));
		self.af_cache = {};
	
	def retrieve_variant(self,xchr,xpos):
		rows_vcf = self.tabix_vcf.fetch(xchr,xpos-1,xpos);
		out_rows = [];
		for row in rows_vcf:
			dict_row = self.row_to_dict(row);
			if int(dict_row['POS']) == xpos:
				out_rows.append(dict_row);
		return(out_rows);

	def row_to_dict(self, row):
		dict_out = {};
		row = row.rstrip().split("\t");
		for xcol, xname in zip(self.dict_cols.keys(),self.dict_cols.values()):
			dict_out[xname] = row[xcol];
		return(dict_out);

	def retrieve_variant_af(self,xchr,xpos,xref,xalt,af_field="AF",cache=False):
		unique_id = "_".join(map(str,[xchr,xpos,xref,xalt]));
		if cache == True and unique_id in self.af_cache:
			return(self.af_cache[unique_id]);
		else:
			rows_vcf = self.tabix_vcf.fetch(xchr,xpos-1,xpos);
			
			for row in rows_vcf:
				dict_row = self.row_to_dict(row);
				if int(dict_row['POS']) == xpos and dict_row['REF']:
					alt_alleles = dict_row['ALT'].split(",");
					if xalt in alt_alleles:
						info_dict = self.dict_from_vcf_info(dict_row['INFO']);
						if af_field in info_dict:
							afs = [];
							for xaf in info_dict[af_field].split(","):
								try:
									afs.append(float(xaf));
								except:
									afs.append(float('nan'));
							alt_index = alt_alleles.index(xalt);
							if alt_index < len(afs):
								af = afs[alt_index];
								if af != float('nan'):
									result = [af,min([1-af,af])];
									if cache == True: self.af_cache[unique_id] = result;
									return(result);

		result = [float('nan'),float('nan')];
		if cache == True: af_cache[unique_id] = result;
		return(result);

	def calculate_variant_af(self,dict_row,xalt,af_field="AF",cache=False):
		unique_id = "_".join(map(str,[dict_row['CHROM'],dict_row['POS'],dict_row['REF'],dict_row['ALT']]));
		if cache == True and unique_id in self.af_cache:
			return(self.af_cache[unique_id]);
		else:
			info_dict = self.dict_from_vcf_info(dict_row['INFO']);
			alt_alleles = dict_row['ALT'].split(",");
			if af_field in info_dict:
				afs = [];
				for xaf in info_dict[af_field].split(","):
					try:
						afs.append(float(xaf));
					except:
						afs.append(float('nan'));
				alt_index = alt_alleles.index(xalt);
				if alt_index < len(afs):
					af = afs[alt_index];
					if af != float('nan'):
						result = [af,min([1-af,af])];
						if cache == True: self.af_cache[unique_id] = result;
						return(result);

	def dict_from_vcf_info(self, info):
		info = info.split(";");
		dict_info = {};
		for item in info:
			fields = item.split("=")
			if len(fields) == 2:
				tag = fields[0];
				value = fields[1];
				dict_info[tag] = value;
		return dict_info;

	def sample_names(self):
		out_names = [];
		for i in range(9,len(self.dict_cols)):
			out_names.append(self.dict_cols[i]);
		return(out_names);

if __name__ == "__main__":
	main();
