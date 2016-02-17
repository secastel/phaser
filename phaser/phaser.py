import multiprocessing;
import string;
import argparse;
import gzip;
import tempfile;
import subprocess;
import itertools;
from intervaltree import IntervalTree;
import sys;
import time;
from scipy.stats import binom;
import numpy;
import os;
import vcf;
import math;

## IMPORTANT NOTES
# ASSUMES THAT THERE ARE 2 REAL HAPLOTYPE COMBINATIONS (DIPLOID, NO TRANS-SPLICING)

## EXTERNAL DEPENDENCIES
# SAMTOOLS
# BEDTOOLS

## PYTHON LIBRARIES
# INTERVALTREE
# PYVCF
# NUMPY
# SCIPY

def main():
	#Arguments passed 
	parser = argparse.ArgumentParser()
	# required
	parser.add_argument("--bam", help="Indexed BAMs (comma separated) containing RNA-seq reads", required = True)
	parser.add_argument("--vcf", help="VCF for sample", required = True)
	parser.add_argument("--sample", default="", help="Sample name in VCF", required = True)
	parser.add_argument("--mapq", help="Minimum MAPQ for reads to be used for phasing. Can be a comma separated list, each value corresponding to the min MAPQ for a file in the input BAM list. Useful in cases when using both for example DNA and RNA libraries which having different mapping qualities.", required = True)
	parser.add_argument("--baseq", type=int, default=0, help="Minimum baseq for bases to be used for phasing", required = True)
	parser.add_argument("--o", help="Out prefix")
	
	# optional
	parser.add_argument("--cc_threshold", type=float, default=0.01, help="Threshold for significant conflicting variant configuration. The connection between any two variants with a conflicting configuration p-value lower than this threshold will be removed.")
	parser.add_argument("--isize", default="0", help="Maximum allowed insert size for read pairs. Can be a comma separated list, each value corresponding to a max isize for a file in the input BAM list. Set to 0 for no maximum size.")
	parser.add_argument("--as_q_cutoff", type=float, default=0.05, help="Bottom quantile to cutoff for read alignment score.")
	parser.add_argument("--ab_q_cutoff", type=float, default=0, help="Bottom quantile to cutoff for read aligned bases.")
	parser.add_argument("--blacklist", default="", help="BED file containing genomic intervals to be excluded from phasing (for example HLA).")
	parser.add_argument("--write_vcf", type=int, default=1, help="Create a VCF containing phasing information (0,1).")
	parser.add_argument("--include_indels", type=int, default=0, help="Include indels in the analysis (0,1). NOTE: since mapping is a problem for indels including them will likely result in poor quality phasing unless specific precautions have been taken.")
	parser.add_argument("--output_read_ids", type=int, default=0, help="Output read IDs in the coverage files (0,1).")
	parser.add_argument("--output_orphans", type=int, default=0, help="Output reads which fail to map to either allele (0, 1).")
	parser.add_argument("--remove_dups", type=int, default=1, help="Remove duplicate reads from all analyses (0,1).")
	parser.add_argument("--pass_only", type=int, default=1, help="Only use variants labled with PASS in the VCF filter field (0,1).")
	parser.add_argument("--min_cov", type=int, default=0, help="Minimum total coverage level before outputting haplotypic counts.")
	parser.add_argument("--unphased_vars", type=int, default=1, help="Output unphased variants (singletons) in the haplotypic_counts and haplotypes files (0,1).")
	
	# genome wide phasing
	parser.add_argument("--gw_phase_method", type=int, default=0, help="Method to use for determing genome wide phasing. NOTE requires input VCF to be phased, and optionally a VCF with allele frequencies (see --gw_af_vcf). 0 = Use most common haplotype phase. 1 = MAF weighted phase anchoring.")
	parser.add_argument("--gw_af_vcf", default="", help="VCF with allele frequencies from the population which was used to do the phasing. If left blank it will look for an allele frequency in the input VCF (--vcf).")
	parser.add_argument("--gw_af_field", default="AF", help="Field from --gw_af_vcf to use for allele frequency.")
	parser.add_argument("--gw_phase_vcf", type=int, default=0, help="Replace GT field of output VCF using phASER genome wide phase (0,1). See --gw_phase_method for options.")
	parser.add_argument("--gw_phase_vcf_min_confidence", type=float, default=0.90, help="If replacing GT field in VCF only replace when phASER haplotype gw_confidence >= this value.")
	
	# performance
	parser.add_argument("--threads", type=int, default=1, help="Maximum number of threads to use. Note the maximum thread count for some tasks is bounded by the data (for example 1 thread per contig for haplotype construction).")
	parser.add_argument("--max_block_size", type=int, default=15, help="Maximum number of variants to phase at once. Number of haplotypes tested = 2 ^ # variants in block. Blocks larger than this will be split into sub blocks, phased, and then the best scoring sub blocks will be phased with each other.")
	parser.add_argument("--reads_mem", type=int, default=0, help="Store reads that overlap variants in memory (0,1). NOTE: storing reads in memory increases speed, but grealtly increases memory overhead.")
	parser.add_argument("--vcf_mem", type=int, default=1, help="Load the entire VCF into memory for quick processing (0,1).")
	parser.add_argument("--temp_dir", default="", help="Location of temporary directory to use for storing files. If left blank will default to system temp dir. NOTE: potentially large files will be stored in this directory, so please ensure there is sufficient free space.")
	parser.add_argument("--max_items_per_thread", type=int, default=100000, help="Maximum number of items that can be assigned to a single thread to process. NOTE: if this number is too high Python will stall when trying to join the pools.")
	
	# debug / development / reporting
	parser.add_argument("--show_warning", type=int, default=0, help="Show warnings in stdout (0,1).")
	parser.add_argument("--debug", type=int, default=0, help="Show debug mode messages (0,1).")
	parser.add_argument("--chr", default="", help="Restrict haplotype phasing to a specific chromosome.")
	parser.add_argument("--unique_ids", type=int, default=0, help="Generate and output unique IDs instead of those provided in the VCF (0,1). NOTE: this should be used if your VCF does not contain a unique ID for each variant.")
	parser.add_argument("--output_network", default="", help="Output the haplotype connection network for the given variant.")
	parser.add_argument("--only_phased", type=int, default=0, help="Only use variants in the input VCF which were marked as phased (0,1).")
	
	global args;
	args = parser.parse_args()
	
	#setup
	version = "0.2";
	fun_flush_print("");
	fun_flush_print("##################################################")
	fun_flush_print("              Welcome to phASER v%s"%(version));
	fun_flush_print("  Author: Stephane Castel (scastel@nygenome.org)")
	fun_flush_print("##################################################");
	fun_flush_print("");
	
	if args.temp_dir != "":
		tempfile.tempdir = args.temp_dir;
	devnull = open(os.devnull, 'w')
	
	# check that all passed files actually exist
	check_files = args.bam.split(",") + [args.vcf,args.blacklist];
	
	for xfile in check_files:
		if xfile != "":
			if os.path.isfile(xfile) == False:
				fatal_error("File: %s not found."%(xfile));
	
	start_time = time.time();
	global dict_blacklist_interval;
	dict_blacklist_interval = {};
	
	if len(args.blacklist) > 0:
		fun_flush_print("#0. Loading blacklist intervals...");
		stream_in = open(args.blacklist,"r");
		for line in stream_in:
			columns = line.replace("\n","").split("\t");
			chr = columns[0];
			if chr not in dict_blacklist_interval: dict_blacklist_interval[chr] = IntervalTree();
			dict_blacklist_interval[chr][int(columns[1]):int(columns[2])] = columns[3];
	
		stream_in.close();
	
	# load the allele frequency VCF if specified
	if args.gw_af_vcf != "":
		if os.path.isfile(args.gw_af_vcf) == True:
			vcf_af = vcf.Reader(filename=args.gw_af_vcf);
		else:
			fatal_error("Allele frequency VCF (--gw_af_vcf) specified does not exit.");
		
	
	fun_flush_print("#1. Loading heterozygous variants into intervals...");
	
	#load VCF
	if ".gz" in args.vcf:
		if args.vcf_mem == 1:
			fun_flush_print("     loading VCF into memory...");
			vcf_in = subprocess.check_output("gunzip -c "+args.vcf, shell=True);
			vcf_lines = vcf_in.split("\n")
		else:
			vcf_in = gzip.open(args.vcf, "r");
	else:
		if args.vcf_mem == 1:
			fun_flush_print("     loading VCF into memory...");
			vcf_in = subprocess.check_output("cat "+args.vcf, shell=True);
			vcf_lines = vcf_in.split("\n")
		else:
			vcf_in = open(args.vcf, "r");
	
	# figure out which column is the correct sample
	global sample_column;
	
	if len(args.sample) > 0:
		sample_column = -1;
	
		if args.vcf_mem == 1:
			vcf_loop = vcf_lines;
		else:
			vcf_loop = vcf_in;
	
		for line in vcf_loop:
			if "#CHROM" in line:
				vcf_columns = line.replace("\n","").split("\t");
				for i in range(len(vcf_columns)):
					if vcf_columns[i] == args.sample:
						sample_column = i;
						break;
			if sample_column != -1:
				break;
	
		if sample_column == -1:
			fatal_error("Sample not found in VCF.");
			quit();
	else:
		#default to first sample if none specified
		sample_column = 9;
		last_chrom = "";
	
	# get all heterozygous sites, load into interval
	global dict_vcf_interval;
	dict_vcf_interval = {};
	global dict_variant_reads;
	dict_variant_reads = {};
	
	if args.vcf_mem == 1:
		fun_flush_print("     parsing VCF...");
		pool_input = pool_split(args.threads, vcf_lines);
		pool_output = parallelize(process_vcf, pool_input);
		del vcf_in;
		del vcf_lines;
		del pool_input;
	else:
		fun_flush_print("     parsing VCF...");
		pool_output = process_vcf(vcf_in)
		vcf_in.close();
	
	bed_out = tempfile.NamedTemporaryFile(delete=False);
	pool_input = {};
	
	total_het_vars = 0;
	total_blacklisted_vars = 0;
	total_indels_excluded = 0;
	
	for item in pool_output:
		total_het_vars += item[0];
		total_blacklisted_vars += item[1];
		total_indels_excluded += item[2];
		out_list = item[3];
		
		for output in out_list:
			# dict_out = {"id":id,"ref":vcf_columns[3],"chr":chrom,"pos":pos,"alleles":ind_alleles,"phase":phase,"maf":maf,"reads":[[] for i in range(len(ind_alleles))]};
			if output != None:
				if output['chr'] not in pool_input: pool_input[output['chr']] = [output['chr']];
		
				dict_variant_reads[output['id']] = output;
				pool_input[output['chr']].append([output['id'],output['pos'],output['length']]);
				#beds are 0 based
				bed_out.write("\t".join([output['chr'],str(output['pos']-1),str(output['pos']),"\n"]));
	
	fun_flush_print("          %d total heterozygous variants, %d indels excluded, %d blacklisted variants"%(total_het_vars,total_indels_excluded,total_blacklisted_vars));
	
	fun_flush_print("     creating genomic intervals...");
	bed_out.close();
	
	pool_output = parallelize(generate_intervals, pool_input.values());
	del pool_input;
	
	for output in pool_output:
		chr = output[0];
		interval = output[1];
		dict_vcf_interval[chr] = interval;
	
	fun_flush_print("#2. Retrieving reads that overlap heterozygous sites...");
	
	#works with multiple input bams
	bam_list = args.bam.split(",");
	
	#mapq
	mapq_list = args.mapq.split(",");
	if len(mapq_list) == 1 and len(bam_list) > 1:
		mapq_list = mapq_list * len(bam_list);
	elif len(mapq_list) != len(bam_list):
		fatal_error("Number of mapq values and input BAMs does not match. Supply either one mapq to be used for all BAMs or one mapq per input BAM.");
	
	#isize
	isize_list = args.isize.split(",");
	if len(isize_list) == 1 and len(bam_list) > 1:
		isize_list = isize_list * len(bam_list);
	elif len(mapq_list) != len(isize_list):
		fatal_error("Number of isize values and input BAMs does not match. Supply either one isize to be used for all BAMs or one isize per input BAM.");
		
	#now get bam reads that overlap het sites using SAMTOOLS
	samtools_arg = "";
	# remove dups if necessary, and only include properly paired read (ie in correct orientation)
	if args.remove_dups == 1:
		samtools_arg = "-F 0x400 -f 2"
	else:
		samtools_arg = "-f 2"
	
	global use_as_cutoff;
	global as_cutoff;
	global ab_cutoff;
	global isize_cutoff;
	
	intersected_read_output = [];
	
	if args.reads_mem == 1:
		#first get all reads that overlap het sites into memory
		for bam, mapq, isize in zip(bam_list, mapq_list, isize_list):
			fun_flush_print("     file: %s"%(bam));
			fun_flush_print("          minimum mapq: %s"%(mapq));
			reads = subprocess.check_output("samtools view "+samtools_arg+"-L "+bed_out.name+" -q "+mapq+" "+bam, stderr=devnull, shell=True);
			reads = reads.split("\n");
			total_reads = len(reads);
			fun_flush_print("          retrieved %d reads"%(total_reads));
		
			use_as_cutoff, as_cutoff = calculate_as_cutoff(reads);
			if use_as_cutoff == True: fun_flush_print("          using alignment score cutoff of %d"%(as_cutoff));
			ab_cutoff = calculate_ab_cutoff(reads);
			if ab_cutoff > 0: fun_flush_print("          using aligned bases cutoff of %d"%(ab_cutoff));
			
			isize_cutoff = int(isize);
			if isize_cutoff > 0: fun_flush_print("          using insert size cutoff of %d"%(isize_cutoff));
			
			pool_input = pool_split(args.threads, reads);
		
			fun_flush_print("          assigning reads to variants...");
			# because reads are in memory can use multiprocessing
			intersected_read_output += parallelize(intersect_reads_vars, pool_input);
			del pool_input;
			del reads;
		
	elif args.reads_mem == 0:
		fun_flush_print("     Reads are being written to disk, this will impact performance.");
		reads_out = tempfile.NamedTemporaryFile(delete=False);
		reads_out.close();
		for bam, mapq, isize in zip(bam_list, mapq_list, isize_list):
			fun_flush_print("     file: %s"%(bam));
			fun_flush_print("          minimum mapq: %s"%(mapq));
			subprocess.call("samtools view "+samtools_arg+" -L "+bed_out.name+" -q "+mapq+" "+bam+" > "+reads_out.name, stderr=devnull, shell=True);
			
			total_reads = int(subprocess.check_output("cat "+reads_out.name+" | wc -l ", stderr=devnull, shell=True));
			fun_flush_print("          retrieved %d reads"%(total_reads));
		
			# read the file in and calculate the AS cutoff
			reads_in = open(reads_out.name, "r");
			use_as_cutoff, as_cutoff = calculate_as_cutoff(reads_in);
			if use_as_cutoff == True: fun_flush_print("          using alignment score cutoff of %d"%(as_cutoff));
			reads_in.close();
		
			# read the file in and calculate the AB cutoff
			reads_in = open(reads_out.name, "r");
			ab_cutoff = calculate_ab_cutoff(reads_in);
			if ab_cutoff > 0: fun_flush_print("          using aligned bases cutoff of %d"%(ab_cutoff));
			reads_in.close();
			
			isize_cutoff = int(isize);
			if isize_cutoff > 0: fun_flush_print("          using insert size cutoff of %d"%(isize_cutoff));
			
			# now split the file for use by multiprocessing
			file_list = [];
			if args.threads > 1:
				data_length = total_reads;
			
				# calculate pool size if all data is divided by number of threads
				optimal_pool_size = data_length/args.threads;
				pool_size = min([args.max_items_per_thread, optimal_pool_size]);
				pool_inputs = data_length / pool_size;
			
				fun_flush_print("          splitting reads into %d files with %d reads"%(pool_inputs, pool_size));
			
				stream_reads_in = open(reads_out.name, "r");
				line_number = 1;
				file_number = 1;
			
				file_name = new_temp_file();
				stream_out = open(file_name, "w");
				file_list.append(file_name);
			
				for line in stream_reads_in:
					if (line_number % pool_size == 0 and file_number != pool_inputs):
						# close current file
						stream_out.close();
					
						# new file
						file_number += 1;
						file_name = new_temp_file();
						stream_out = open(file_name, "w");
						file_list.append(file_name);
					
					# write to open file
					stream_out.write(line);
					line_number += 1;
				
				# close current file
				stream_out.close();
			else:
				file_list.append(reads_out.name);
		
			fun_flush_print("          assigning reads to variants...");
			intersected_read_output += parallelize(intersect_reads_vars_file, file_list);
		
		#cleanup temp files
		os.remove(reads_out.name);
		for file in file_list:
			if os.path.isfile(file):
				os.remove(file);
	
	os.remove(bed_out.name);
	
	# process pool outputs
	read_intersects = [];
	orphan_reads = "";
	base_matches = {};
	base_mismatches = {};
	variant_overlaps = [];
	
	for item in intersected_read_output:
		read_intersects += item[0];
		orphan_reads += item[1];
		base_matches = { k: base_matches.get(k, 0) + item[2].get(k, 0) for k in set(base_matches) | set(item[2]) }
		base_mismatches = { k: base_mismatches.get(k, 0) + item[3].get(k, 0) for k in set(base_mismatches) | set(item[3]) }
		variant_overlaps += item[4];
	
	#clear memory
	del intersected_read_output;
	
	# calculate noise level
	base_match_count = 0;
	base_mismatch_count = 0;
	
	for variant in base_matches:
		mis_matches = 0;
		
		if variant in base_mismatches:
			mis_matches = base_mismatches[variant];
		
		# require other bases to be < 5% of total coverage for this variant
		# protects against genotyping errors
		if base_matches[variant] > 0 and (float(mis_matches) / float(mis_matches+base_matches[variant])) < 0.50:
			base_match_count += base_matches[variant];
			base_mismatch_count += mis_matches;
	
	if base_match_count == 0:
		fatal_error("No reads could be matched to variants. Please double check your settings and input files. Common reasons for this occurring include: 1) MAPQ or BASEQ set incorrectly 2) BAM and VCF have different chromosome names (IE 'chr1' vs '1').");
	
	fun_flush_print("#3. Identifying connected variants...");
	# probability of generating a random base
	global noise_e;
	noise_e = (float(base_mismatch_count) / (float(base_match_count+base_mismatch_count)*2));
	fun_flush_print("     sequencing noise level estimated at %f"%(noise_e));
	
	# now add the reads to the variant dictionaries
	read_vars = {};
	for intersect in read_intersects:
		var_id = intersect[0];
		allele_index = intersect[1];
		read_id = intersect[2];
		var_dict = dict_variant_reads[var_id];
		if allele_index >= 0:
			var_dict['reads'][allele_index].append(read_id);
		else:
			var_dict['other_reads'].append(read_id);
		if read_id not in read_vars: read_vars[read_id] = [];
		if var_id not in read_vars[read_id]: read_vars[read_id].append(var_id);
	
	# premake read sets for faster comparison
	for var_id in dict_variant_reads:
		dict_variant_reads[var_id]['read_set'] = [];
		for allele_reads in dict_variant_reads[var_id]['reads']:
			dict_variant_reads[var_id]['read_set'].append(set(allele_reads));
		dict_variant_reads[var_id]['other_read_set'] = set(dict_variant_reads[var_id]['other_reads']);

	# now create the quick lookup dictionary
	# this is used for haplotype construction
	# dictionary tells you what variants are connected
	dict_variant_overlap = {};
	for read_id in read_vars.keys():
		overlapped_variants = read_vars[read_id];
		for variant in overlapped_variants:
			var_chr = dict_variant_reads[variant]['chr'];
			for other_variant in overlapped_variants:
				other_var_chr = dict_variant_reads[other_variant]['chr'];
				# Restrict to being on the same chromosome, speeds up and allows parallelization
				# might not be desired for some very specific cases (ie trans-splicing)
				if var_chr == other_var_chr and other_variant != variant:
					if var_chr not in dict_variant_overlap: dict_variant_overlap[var_chr] = {};
					if variant not in dict_variant_overlap[var_chr]: dict_variant_overlap[var_chr][variant] = [];
					dict_variant_overlap[var_chr][variant].append(other_variant);
	
	# make sets of overlaps
	for chr in dict_variant_overlap:
		for variant in dict_variant_overlap[chr]:
			dict_variant_overlap[chr][variant] = set(dict_variant_overlap[chr][variant]);
	
	## now run the test to determine if the number of reads with conflicting connections is
	## higher than noise for a given variant pair.
	## if so these two variants will be disconnected, so that they won't be used for haplotype construction
	tested_connections = set([]);
	pool_input = [];
	
	for chr in dict_variant_overlap:
		for variant_a in dict_variant_overlap[chr]:
			overlapping_variants = dict_variant_overlap[chr][variant_a];
			for variant_b in overlapping_variants:
				key1 = variant_a+"|"+variant_b;
				key2 = variant_b+"|"+variant_a;
				if key1 not in tested_connections and key2 not in tested_connections:
					pool_input.append([chr,variant_a,variant_b]);
					tested_connections.add(key1);
	
	pool_output = parallelize(test_variant_connection, pool_input);
	
	out_stream = open(args.o+".variant_connections.txt","w");
	out_stream.write("variant_a\tvariant_b\tsupporting_connections\ttotal_connections\tconflicting_configuration_p\tphase_concordant\n");
	
	# remove all those connections which failed
	c_dropped = 0;
	for connection in pool_output:
		chr,variant_a,variant_b,conflicting_config_p,c_supporting,c_total,phase_concordant = connection;
		
		# if the number of conflicting reads is more than would be expected from noise, then disconnect these two variants
		# they will not be used for haplotype construction
		out_stream.write("\t".join(map(str,[variant_a,variant_b,c_supporting,c_total,conflicting_config_p,phase_concordant]))+"\n");
		if conflicting_config_p < args.cc_threshold:
			#print("%s	%s"%(variant_a,variant_b));
			dict_variant_overlap[chr][variant_a].remove(variant_b);
			dict_variant_overlap[chr][variant_b].remove(variant_a);
			
			# if these variants have no other connections remove them from overlap dictionary
			if len(dict_variant_overlap[chr][variant_a]) == 0:
				del dict_variant_overlap[chr][variant_a];
			if len(dict_variant_overlap[chr][variant_b]) == 0:
				del dict_variant_overlap[chr][variant_b];
			
			c_dropped += 1;
	
	out_stream.close();
	
	fun_flush_print("     %d variant connections dropped because of conflicting configurations (threshold = %f)"%(c_dropped,args.cc_threshold));
		
	# create output file for orphan reads if needed
	if args.output_orphans == 1:
		orphan_out = open(args.o+".orphans.sam","w");
		orphan_out.write(orphan_reads);
		orphan_out.close();
	
	#clear memory
	del orphan_reads;
	
	# output the coverage level per snp
	# same format as GATK tool:
	stream_out = open(args.o+".allelic_counts.txt","w");
	stream_out.write("contig	position	variantID	refAllele	altAllele	refCount	altCount	totalCount\n");
	covered_count = 0;
	
	for variant in dict_variant_reads:
		snp_dict = dict_variant_reads[variant];
		
		ref_reads = len(set(snp_dict['reads'][0]));
		alt_reads = len(set(snp_dict['reads'][1]));
		if ref_reads+alt_reads > 0:
			covered_count += 1;
			stream_out.write("\t".join([snp_dict['chr'],str(snp_dict['pos']),variant,snp_dict['alleles'][0],snp_dict['alleles'][1],str(ref_reads),str(alt_reads),str(ref_reads+alt_reads)+"\n"]));
	stream_out.close();
	
	fun_flush_print("     %d variants covered by at least 1 read"%(covered_count));
	
	# close disk space storage out
	if args.reads_mem == 0:
		reads_out.close();	
	
	# record the total number of hets
	total_het_variants = len(dict_variant_reads.keys());
	
	remove_keys = [];
	if args.unphased_vars == 0:
		# clear all SNPs with no connections to others from the dictionary to free up memory
		# if we don;t want to output unphased snps
		
		for variant in dict_variant_reads:
			chr = dict_variant_reads[variant]['chr'];
			if chr not in dict_variant_overlap:
				remove_keys.append(variant);
			elif variant not in dict_variant_overlap[chr]:
				remove_keys.append(variant);
	else:
		# otherwise just remove variants with 0 coverage
		for variant in dict_variant_reads:
			if len(dict_variant_reads[variant]['reads'][0]) + len(dict_variant_reads[variant]['reads'][1]) == 0:
				remove_keys.append(variant);
	
	for key in remove_keys:
		del dict_variant_reads[key];
		
	print_debug("     removed %d variants from memory in cleanup"%(len(remove_keys)));
	
	# using only the overlapping SNP dictionary build haplotype blocks
	fun_flush_print("#4. Identifying haplotype blocks...");
	
	block_haplotypes = [];
	phased_vars = 0;
	
	pool_output = parallelize(build_haplotypes, dict_variant_overlap.values());
	
	for chr_haplotypes in pool_output:
		for haplotype_block in chr_haplotypes:
			block_haplotypes.append(haplotype_block);
			phased_vars += len(haplotype_block);
	
		
	# count total number of junctions for each haplotype
	pool_output = parallelize(count_hap_junctions, block_haplotypes);
	
	haplotype_junction_counts = {};
	for output in pool_output:
		key = ",".join(output[0]);
		haplotype_junction_counts[key] = output[1];
	
	# now for each of the blocks identify the phasing with the most supporting reads
	fun_flush_print("#5. Phasing blocks...");
	pool_input = [];
	large_blocks = [];
	
	i = 1;
	for block in block_haplotypes:
		# produce all possible allele configurations for length X variants
		i += 1;
		if len(block) <= args.max_block_size:
			configurations = ["".join(seq) for seq in itertools.product("01", repeat=len(block))];
			for configuration in configurations:
				pool_input.append([block,configuration]);
		else:
			# process these large blocks separately
			large_blocks.append(block);
			print_warning("     maximum block size exceeded for block: "+str(block));
		
	pool_output = parallelize(count_hap_reads, pool_input);
	del pool_input;
	
	supporting_reads = {};
	for item in pool_output:
		block = list_to_string(item[0],",");
		configuration = item[1];
		count = item[2];
		if block not in supporting_reads: supporting_reads[block] = {};
		
		supporting_reads[block][configuration] = count;
	
	# run the haplotype construction for large blocks
	pool_input = [];
	if len(large_blocks) > 0:
		fun_flush_print("     phasing large (>%d variants) blocks..."%(args.max_block_size));
		
		# break the large blocks into smaller chunks (of max block size)
		# then phase the smaller chunks with each other
		# variants must be sorted by position
		
		for block in large_blocks:
			split_blocks = [];
	
			n_sub_blocks = len(block) / args.max_block_size;
			remainder = len(block) % args.max_block_size;
			for i in range(0,n_sub_blocks):
				split_blocks.append(block[i*args.max_block_size:(i+1)*args.max_block_size]);
	
			#remainder
			if remainder > 0:
				split_blocks.append(block[(i+1)*args.max_block_size:]);
	
			blockn = 0;
			for sblock in split_blocks:
				configurations = ["".join(seq) for seq in itertools.product("01", repeat=len(sblock))];
				for configuration in configurations:
					pool_input.append([sblock,configuration,block,blockn]);
				blockn += 1;
					
		# count haplotype reads in each sub block
		pool_output = parallelize(count_hap_reads, pool_input);
		del pool_input;
		
		# now for each sub block find the most supported haplotype
		dict_parent_blocks = {};
		
		for item in pool_output:
			sblock = list_to_string(item[0],",");
			configuration = item[1];
			count = item[2];
			pblock = list_to_string(item[3]);
			blockn = item[4];
			
			if pblock not in dict_parent_blocks: dict_parent_blocks[pblock] = {};
			if sblock not in dict_parent_blocks[pblock]: dict_parent_blocks[pblock][sblock] = {'BLOCK_N':blockn};
		
			dict_parent_blocks[pblock][sblock][configuration] = count;
		
		pool_input = [];
		for parent_block in dict_parent_blocks:
			n_sub_blocks = len(block) / args.max_block_size;
			remainder = len(block) % args.max_block_size
			if remainder > 0: n_sub_blocks += 1;
			
			best_configs = {};
			
			for sub_block in dict_parent_blocks[parent_block]:
				blockn = dict_parent_blocks[parent_block][sub_block]['BLOCK_N'];
				key = str(blockn)+","+sub_block;
				best_configs[key] = [];
				
				configurations = set(dict_parent_blocks[parent_block][sub_block].keys()) - set(['BLOCK_N']);
				resolved_configs = {};
				while len(configurations) > 0:
					hap_a = configurations.pop();
					hap_b = "".join([str(int(not bool(int(x)))) for x in hap_a]);
					resolved_configs[hap_a+"|"+hap_b] = dict_parent_blocks[parent_block][sub_block][hap_a] + dict_parent_blocks[parent_block][sub_block][hap_b];
					configurations.remove(hap_b);
		
				# take all configs tying for max count
				max_count = max(resolved_configs.values());
				for config in resolved_configs:
					if resolved_configs[config] == max_count:
						best_configs[key] += config.split("|");
			
			# sort by block order
			xsplit = [x.split(",") for x in best_configs.keys()];
			xsort = sorted(xsplit, key = lambda x: (int(x[0])));
			
			# now generate only combinations of the best possible sub_blocks to be tested
			haplo_vars = [];
			haplo_configs = [];
			
			for vars in xsort:
				key = ",".join(vars);
				haplo_vars += vars[1:];
				haplo_configs.append(best_configs[key]);
				
			test_configs = ["".join(x) for x in list(itertools.product(*haplo_configs))];
			
			for config in test_configs:
				pool_input.append([haplo_vars,config]);
		
		# now phase the best blocks
		pool_output = parallelize(count_hap_reads, pool_input);
		del pool_input;
		
		for item in pool_output:
			block = list_to_string(item[0],",");
			configuration = item[1];
			count = item[2];
			if block not in supporting_reads: supporting_reads[block] = {};
		
			supporting_reads[block][configuration] = count;
	
	
	# now determine the haplotype with the most support
	fun_flush_print("     identifying haplotypes with most support...");
	final_haplotypes = [];
	for block in supporting_reads:
		configurations = set(supporting_reads[block].keys());
		resolved_configs = {};
		while len(configurations) > 0:
			hap_a = configurations.pop();
			hap_b = "".join([str(int(not bool(int(x)))) for x in hap_a]);
			resolved_configs[hap_a+"|"+hap_b] = supporting_reads[block][hap_a] + supporting_reads[block][hap_b];
			configurations.remove(hap_b);
		
		max_count = max(resolved_configs.values());
		if resolved_configs.values().count(max_count) > 1:
			print_warning("     inconclusive phasing for %s"%(block));
		else:
			for config in resolved_configs:
				if resolved_configs[config] == max_count:
					#print("     most supported phase for %s = %s"%(block, config));
					final_haplotypes.append([block.split(","),config,resolved_configs[config]])
		
	fun_flush_print("#6. Outputting haplotypes...");	
	
	stream_out_ase = open(args.o+".haplotypic_counts.txt","w");
	ase_columns = ["contig","start","stop","variants","variantCount","haplotypeA","haplotypeB","aCount","bCount","totalCount","blockGWPhase","gwStat"];
	if args.output_read_ids == 1:
		ase_columns += ["read_ids_a","read_ids_b"];
	stream_out_ase.write("\t".join(ase_columns)+"\n");
	
	stream_out = open(args.o+".haplotypes.txt","w");
	stream_out.write("\t".join(['contig','start','stop','length','variants','variant_ids','variant_alleles','reads_hap_a','reads_hap_b','reads_total','connections_supporting','connections_total','annotated_phase','phase_concordant','gw_phase','gw_confidence'])+"\n");
	
	stream_out_allele_configs = open(args.o+".allele_config.txt","w");
	stream_out_allele_configs.write("\t".join(['variant_a','rsid_a','variant_b','rsid_b','configuration'])+"\n");
		
	global haplotype_lookup;
	haplotype_lookup = {};
	global haplotype_pvalue_lookup;
	haplotype_pvalue_lookup = {};
	global haplotype_gw_stat_lookup;
	haplotype_gw_stat_lookup = {};
	
	all_variants = [];
	
	block_index = 0;
	
	for block in final_haplotypes:
		#get all unique variants
		block_index += 1;
		variants = block[0];
		all_variants += variants;
		
		haplotype_a = block[1].split("|")[0];
		haplotype_b = block[1].split("|")[1];
		supporting_connections = block[2];
		total_connections = haplotype_junction_counts[list_to_string(variants)];
		
		if args.unique_ids == 0:
			rsids = [dict_variant_reads[x]['rsid'] for x in variants];
		else:
			rsids = variants;
				
		chrs = [dict_variant_reads[x]['chr'] for x in variants];
		positions = map(int, [dict_variant_reads[x]['pos'] for x in variants]);
		
		# if at some point we want to implement a haplotype test do it here
		hap_p = 0;
		haplotype_pvalue_lookup[list_to_string(variants)] = hap_p;
		
		for var_index in range(0, len(variants)):	
			id = variants[var_index];
			haplotype_lookup[id] = [variants, haplotype_a[var_index]+"|"+haplotype_b[var_index],block_index];
			
		alleles = [[],[]];
		phases = [[],[]];
		set_reads = [[],[]];
		hap_counts = [0,0];
		
		for hap_index in range(0,2):
			hap_x = [haplotype_a, haplotype_b][hap_index];
			
			for var_index in range(0, len(variants)):
				id = variants[var_index];
				allele = dict_variant_reads[id]['alleles'][int(hap_x[var_index])];
				alleles[hap_index].append(allele);
				phases[hap_index].append(get_allele_phase(allele,dict_variant_reads[id]));
				
				allele_index = dict_variant_reads[id]['alleles'].index(allele);
				set_reads[hap_index] += dict_variant_reads[id]['reads'][allele_index];
				
			set_reads[hap_index] = list(set(set_reads[hap_index]));
			hap_counts[hap_index] = len(set_reads[hap_index]);
			
		# determine if phasing is completely concordant
		# don't include variants whose phase was unknown in the original VCF
		use_phases = filter(lambda x: str(x) != "nan", phases[0]);
		if len(set(use_phases)) <= 1:
			phase_concordant = 1;
		else:
			phase_concordant = 0;
		
		phase_string = ["",""]
		phase_string[0] = "".join([str(x).replace("nan", "-") for x in phases[0]]);
		phase_string[1] = "".join([str(x).replace("nan", "-") for x in phases[1]]);
		
		### GENOME WIDE PHASING
		# how many population phased variants do we have in this hap
		nan_strip = [int(x) for x in phases[0] if x >= 0];
		
		# by default corrected is the same as population
		corrected_phases = [phases[0],phases[1]];
		cor_phase_stat = 0.5;
		
		if len(nan_strip) > 0:
			# if setting is on determine genome wide phasing
			if args.gw_phase_method == 0:
				# phase using most common phase
					cor_phase_stat = numpy.mean(nan_strip);
			
					if cor_phase_stat < 0.5:
						corrected_phases = [[0]*len(variants),[1]*len(variants)];
					elif cor_phase_stat > 0.5:
						corrected_phases = [[1]*len(variants),[0]*len(variants)];
					else:
						# no consensus, use population phasing
						print_warning("No GW phasing consensus for %s using method 1"%(str(variants)));
					
					cor_phase_stat = max([cor_phase_stat, 1-cor_phase_stat]);
					
			elif args.gw_phase_method == 1:
				# phase using MAF weighted phase
				# we need the mafs for this, so we need to look them up
				# Step 2 get allele frequencies
				# first get the allele frequency for each of the variants
				haplotype_mafs = [];
				for variant in variants:
					if os.path.isfile(args.gw_af_vcf) == True:
						try:
							var_chrom = dict_variant_reads[variant]['chr'];
							var_pos = dict_variant_reads[variant]['pos'];
							records = vcf_af.fetch(var_chrom,var_pos-1,var_pos);
							for record in records:
								if record.POS == var_pos:
									# get MAF for each allele;
									var_alleles = list(map(str,record.ALT));
									var_afs = map(float, record.INFO[args.gw_af_field]);
									var_mafs = [];
									for var_allele, var_af in zip(var_alleles, var_afs):
										if var_allele in dict_variant_reads[variant]['alleles']:
											var_mafs.append(min([1-var_af,var_af]));
							if len(var_mafs) > 0:
								haplotype_mafs.append(min(var_mafs));
							else:
								haplotype_mafs.append(float('nan'));
						except:
							print_warning("AF Lookup failed (2) for %s:%d"%(chr,pos));
					else:
						if dict_variant_reads[variant]['maf'] != None:
							haplotype_mafs.append(dict_variant_reads[variant]['maf']);
						else:
							haplotype_mafs.append(float('nan'));
			
				if len(haplotype_mafs) == len(variants):
					phase_support = [0,0];
					# now we need to weight the phasing by their MAF
					for phase, maf in zip(phases[0],haplotype_mafs):
						if phase == 0:
							phase_support[0] += maf;
						elif phase == 1:
							phase_support[1] += maf;
				
					# now select the phase with the most MAF support
					if phase_support > 0:
						cor_phase_stat = max(phase_support) / sum(phase_support);
					else:
						cor_phase_stat = 0.5;
				
					if phase_support[0] > phase_support[1]:
						corrected_phases = [[0]*len(variants),[1]*len(variants)];
					elif phase_support[1] > phase_support[0]:
						corrected_phases = [[1]*len(variants),[0]*len(variants)];
					else:
						# no consensus, use population phasing
						print_warning("No GW phasing consensus for %s using method 2"%(str(variants)));
				else:
					print_warning("GW phasing failed for %s"%(str(variants)));
				
		# save the stat for lookup when generating VCF
		haplotype_gw_stat_lookup[list_to_string(variants)] = cor_phase_stat;
		
		# update the variants with their corrected phases
		for var_index in range(0,len(variants)):
			variant = variants[var_index];
			allele_index = dict_variant_reads[variant]['alleles'].index(alleles[0][var_index])
			dict_variant_reads[variant]['gw_phase'][allele_index] = corrected_phases[0][var_index];
			dict_variant_reads[variant]['gw_phase'][1-allele_index] = corrected_phases[1][var_index];
		
		corrected_phase_string = ["",""]
		corrected_phase_string[0] = "".join([str(x).replace("nan", "-") for x in corrected_phases[0]]);
		corrected_phase_string[1] = "".join([str(x).replace("nan", "-") for x in corrected_phases[1]]);
		
		## write the haplotype details
		stream_out.write(str_join("\t",[chrs[0],min(positions),max(positions),max(positions)-min(positions),len(variants),list_to_string(rsids),list_to_string(alleles[0])+"|"+list_to_string(alleles[1]),hap_counts[0],hap_counts[1],sum(hap_counts),supporting_connections,total_connections,phase_string[0]+"|"+phase_string[1],phase_concordant,corrected_phase_string[0]+"|"+corrected_phase_string[1],cor_phase_stat])+"\n");
		
		#$ write ASE stats
		total_cov = int(hap_counts[0])+int(hap_counts[1]);
		if 	total_cov >= args.min_cov:
			out_block_gw_phase = "0/1";
			if corrected_phases[0][0] == 0:
				# haplotype A = GW phase 0
				out_block_gw_phase = "0|1";
			elif corrected_phases[0][0] == 1:
				# haplotype A = GW phase 1
				out_block_gw_phase = "1|0";
			
			fields_out = [chrs[0],min(positions),max(positions),list_to_string(variants),len(variants),list_to_string(alleles[0]),list_to_string(alleles[1]),hap_counts[0],hap_counts[1],total_cov,out_block_gw_phase,cor_phase_stat];
			if args.output_read_ids == 1:
				fields_out += [list_to_string(set_reads[0]),list_to_string(set_reads[1])];
			stream_out_ase.write(str_join("\t",fields_out)+"\n");
		
		## OUTPUT THE NETWORK FOR A SPECIFIC HAPLOTYPE
		if args.output_network in variants:
			#hap_a_network = generate_hap_network([variants, haplotype_a])[0];
			#hap_b_network = generate_hap_network([variants, haplotype_b])[0];
			hap_a_network = generate_hap_network_all(variants)[0];
			stream_out_network = open(args.o+".network.links.txt","w");
			stream_out_network.write("\t".join(["variantA","variantB","connections","inferred\n"]));
			nodes = [];
			#for item in hap_a_network + hap_b_network:
			hap_a_vars = [];
			for item in hap_a_network:
				if item[2] > 0:
					stream_out_network.write(list_to_string(item,"\t")+"\n");
					nodes.append(item[0]);
					nodes.append(item[1]);
			
			stream_out_network.close();
			stream_out_network = open(args.o+".network.nodes.txt","w");
			stream_out_network.write("id\tindex\tassigned_hap\n");
			for item in set(nodes):
				xvar = item.split(":")[0];
				xallele = item.split(":")[1];
				var_index = variants.index(xvar);
				if alleles[0][var_index] == xallele:
					assigned_hap = "A";
				else:
					assigned_hap = "B";
				stream_out_network.write(item+"\t"+str(var_index)+"\t"+assigned_hap+"\n");
			stream_out_network.close()
		
		## OUTPUT allele configuration
		for variant_a, allele_a in zip(variants, alleles[0]):
			for variant_b, allele_b in zip(variants, alleles[1]):
				if variant_a != variant_b:
					a_config = "";
					if (dict_variant_reads[variant_a]['ref'] == allele_a and dict_variant_reads[variant_b]['ref'] == allele_b) or (dict_variant_reads[variant_a]['ref'] != allele_a and dict_variant_reads[variant_b]['ref'] != allele_b):
						# ref and ref are in trans
						# this is a compound het
						a_config = "trans";
					elif (dict_variant_reads[variant_a]['ref'] == allele_a and dict_variant_reads[variant_b]['ref'] != allele_b) or (dict_variant_reads[variant_a]['ref'] != allele_a and dict_variant_reads[variant_b]['ref'] == allele_b):
						# ref and ref are in cis
						a_config = "cis";
					if a_config != "":
						stream_out_allele_configs.write("\t".join([variant_a,dict_variant_reads[variant_a]['rsid'],variant_b,dict_variant_reads[variant_b]['rsid'],a_config])+"\n");
		
	#output read counts for unphased variants
	if args.unphased_vars == 1:
		singletons = set(dict_variant_reads.keys()) - set(all_variants);
	
		for variant in singletons:
			dict_var = dict_variant_reads[variant];
			total_cov = len(dict_var['read_set'][0])+len(dict_var['read_set'][1]);
			if 	total_cov >= args.min_cov:
				if "-" not in dict_var['phase']:
					phase_string = str(dict_var['phase'].index(dict_var['alleles'][0]))+"|"+str(dict_var['phase'].index(dict_var['alleles'][1]));
				else:
					phase_string = "0/1";
				fields_out = [dict_var['chr'],str(dict_var['pos']),str(dict_var['pos']),variant,str(1),dict_var['alleles'][0],dict_var['alleles'][1],str(len(dict_var['read_set'][0])),str(len(dict_var['read_set'][1])),str(total_cov),phase_string,"1"];
				
				if args.output_read_ids == 1:
					fields_out += [list_to_string(dict_var['read_set'][0]),list_to_string(dict_var['read_set'][1])];
			
				stream_out_ase.write("\t".join(fields_out)+"\n");
	
		#output haplotypes for unphased variants (if enabled)
		for variant in singletons:
			dict_var = dict_variant_reads[variant];
			total_cov = len(dict_var['read_set'][0])+len(dict_var['read_set'][1]);
			
			# make sure it is actually phased
			if "-" not in dict_var['phase']:
				phase_string = str(dict_var['phase'].index(dict_var['alleles'][0]))+"|"+str(dict_var['phase'].index(dict_var['alleles'][1]));
			else:
				phase_string = "-|-";
			
			stream_out.write(dict_var['chr']+"\t"+str(dict_var['pos']-1)+"\t"+str(dict_var['pos'])+"\t"+str(1)+"\t"+str(1)+"\t"+variant+"\t"+dict_var['alleles'][0]+"|"+dict_var['alleles'][1]+"\t"+str(len(dict_var['read_set'][0]))+"\t"+str(len(dict_var['read_set'][1]))+"\t"+str(total_cov)+"\t"+str(0)+"\t"+str(0)+"\t"+phase_string+"\t"+str(float('nan'))+"\t"+phase_string+"\t"+str(float('nan'))+"\n");
			
	stream_out.close();
	stream_out_ase.close();
	stream_out_allele_configs.close();
	
	# output VCF
	if args.write_vcf == 1:
		write_vcf();
		
	total_time = time.time() - start_time;
	
	fun_flush_print("COMPLETED - phased %d of %d variants (= %f) using %d reads in %d seconds using %d threads."%(len(all_variants),total_het_variants,float(len(all_variants))/float(total_het_variants),total_reads,total_time,args.threads));

def test_variant_connection(input):
	global noise_e;
	global dict_variant_reads;
	
	chr, variant_a, variant_b = input;
	
	# there are only two possible configurations, determine evidence for each
	# a[ref]b[ref] | a[alt]b[alt]
	hap_config_a_support = len(dict_variant_reads[variant_a]['read_set'][0] & dict_variant_reads[variant_b]['read_set'][0]) + len(dict_variant_reads[variant_a]['read_set'][1] & dict_variant_reads[variant_b]['read_set'][1])
	# a[ref]b[alt] | a[alt]b[ref]
	hap_config_b_support = len(dict_variant_reads[variant_a]['read_set'][1] & dict_variant_reads[variant_b]['read_set'][0]) + len(dict_variant_reads[variant_a]['read_set'][0] & dict_variant_reads[variant_b]['read_set'][1])
		
	# determine if phasing is concordant with what as specified in the input VCF
	phase_concordant = ".";
	# make sure the input VCF had phase
	if "-" not in dict_variant_reads[variant_a]['phase'] and "-" not in dict_variant_reads[variant_b]['phase']:
		if hap_config_a_support > hap_config_b_support:
			
			if dict_variant_reads[variant_a]['phase'].index(dict_variant_reads[variant_a]['alleles'][0]) == dict_variant_reads[variant_b]['phase'].index(dict_variant_reads[variant_b]['alleles'][0]):
				phase_concordant = 1;
			else:
				phase_concordant = 0;
		elif hap_config_a_support < hap_config_b_support:
			if dict_variant_reads[variant_a]['phase'].index(dict_variant_reads[variant_a]['alleles'][1]) == dict_variant_reads[variant_b]['phase'].index(dict_variant_reads[variant_b]['alleles'][0]):
				phase_concordant = 1;
			else:
				phase_concordant = 0;

	# also get the connections from reads where the bases did not match to either ref or alt
	# a[other] -> b[ref]
	other_base_connections = len(dict_variant_reads[variant_a]['other_read_set'] & dict_variant_reads[variant_b]['read_set'][0]);
	# a[other] -> b[alt]
	other_base_connections += len(dict_variant_reads[variant_a]['other_read_set'] & dict_variant_reads[variant_b]['read_set'][1]);
	# a[ref] -> b[other]
	other_base_connections += len(dict_variant_reads[variant_a]['read_set'][0] & dict_variant_reads[variant_b]['other_read_set']);
	# a[alt] -> b[other]
	other_base_connections += len(dict_variant_reads[variant_a]['read_set'][1] & dict_variant_reads[variant_b]['other_read_set']);
	# a[other] -> b[other]
	other_base_connections += len(dict_variant_reads[variant_a]['other_read_set'] & dict_variant_reads[variant_b]['other_read_set']);
	
	c_supporting = max(hap_config_a_support,hap_config_b_support);
	c_total = hap_config_a_support + hap_config_b_support + other_base_connections;
	
	# if no reads support the phase then strip this connection
	if c_supporting == 0:
		conflicting_config_p = 0;
	elif c_total - c_supporting > 0:
		# otherwise do the test
		conflicting_config_p = binom.cdf(c_supporting,c_total,1-((6*noise_e)+(10*math.pow(noise_e,2))));
	else:
		# only both doing the test if there are some conflicting reads
		conflicting_config_p = 1;
		
	return([chr,variant_a,variant_b,conflicting_config_p,c_supporting,c_total,phase_concordant]);

def new_temp_file():
	xfile = tempfile.NamedTemporaryFile(delete=False)
	xfile.close();
	return(xfile.name);

def write_vcf():
	global args;
	global haplotype_lookup;
	global dict_variant_reads;
	global haplotype_pvalue_lookup
	
	fun_flush_print("#7. Outputting phased VCF...");	
	if ".gz" in args.vcf:
		vcf_in = gzip.open(args.vcf,"r");
		vcf_out = gzip.open(args.o+".vcf.gz","w");
	else:
		vcf_in = open(args.vcf,"r");
		vcf_out = open(args.o+".vcf","w");
	
	phase_corrections = 0;
	
	set_phased_vars = set(haplotype_lookup.keys());
	format_section = False;
	for line in vcf_in:
		if "##FORMAT" in line:
			format_section = True;
		elif format_section == True:
			# insert new format info here
			vcf_out.write("##FORMAT=<ID=PG,NUMBER=1,TYPE=String,Description=\"phASER Local Genotype\">\n");
			vcf_out.write("##FORMAT=<ID=PB,NUMBER=1,TYPE=String,Description=\"phASER Local Block\">\n");
			vcf_out.write("##FORMAT=<ID=PI,NUMBER=1,TYPE=String,Description=\"phASER Local Block Index (unique for each block)\">\n");
			vcf_out.write("##FORMAT=<ID=PW,NUMBER=1,TYPE=String,Description=\"phASER Genome Wide Genotype\">\n");
			vcf_out.write("##FORMAT=<ID=PC,NUMBER=1,TYPE=String,Description=\"phASER Genome Wide Confidence\">\n");
			
			format_section = False;
			
		vcf_columns = line.replace("\n","").split("\t");
		if line[0:1] == "#":
			vcf_out.write(line);
			if line[0:4] == "#CHR":
				if args.sample in vcf_columns:
					sample_column = vcf_columns.index(args.sample);
		else:
			##CHROM  POS     ID      REF     ALT     QUAL    FILTER  INFO    FORMAT  NA06986
			id = vcf_columns[2];
			chrom = vcf_columns[0];
			pos = int(vcf_columns[1]);
	
			if "GT" in vcf_columns[8]:
				gt_index = vcf_columns[8].split(":").index("GT");
				genotype = list(vcf_columns[sample_column].split(":")[gt_index]);
				
				if "|" in genotype: genotype.remove("|");
				if "/" in genotype: genotype.remove("/");
			
				# get only the alleles this individual has
				alt_alleles = vcf_columns[4].split(",");
				all_alleles = [vcf_columns[3]] + alt_alleles;
				ind_alleles = [];

				for i in range(0,len(all_alleles)):
					if str(i) in genotype:
						ind_alleles.append(all_alleles[i]);
				
				# make sure there are as many entries in each sample as there should be before adding new columns
				# if there are entries missing add blanks
				n_fields = len(vcf_columns[8].split(":"));
				for i in range(9, len(vcf_columns)):
					sample_fields = len(vcf_columns[i].split(":"));
					if sample_fields != n_fields:
						missing_cols = n_fields - sample_fields;
						vcf_columns[i] += ":" * missing_cols;
				
				vcf_columns[8] += ":PG:PB:PI:PW:PC";
			
				#generate a unique id
				unique_id = chrom+"_"+str(pos)+"_"+("_".join(all_alleles));
					
				if unique_id in set_phased_vars:					
					# retrieve the correct allele number of each allele
					# issue because if a site is multi-allelic it will be converted to 0/1 (ie 0/2 converted to 0/1)
					alleles_out = [];
					gw_phase_out = ["",""];
					block_index = haplotype_lookup[unique_id][2];
					
					for allele in haplotype_lookup[unique_id][1].split("|"):
						allele_base = dict_variant_reads[unique_id]['alleles'][int(allele)];
						vcf_allele_index = all_alleles.index(allele_base);
						
						# get the genome wide phase
						gw_phase = dict_variant_reads[unique_id]['gw_phase'][int(allele)]
						if isinstance(gw_phase, int) == True:
							gw_phase_out[gw_phase] = str(vcf_allele_index);
							
						alleles_out.append(str(vcf_allele_index));
					
					# get rsID for each of the variants on the haplotype
					variants_out = [];
					for variant in haplotype_lookup[unique_id][0]:
						variants_out.append(dict_variant_reads[variant]['rsid']);
					# get the p-value, if there was one for the block
					# pval = haplotype_pvalue_lookup[list_to_string(haplotype_lookup[unique_id][0])];
					gw_stat = haplotype_gw_stat_lookup[list_to_string(haplotype_lookup[unique_id][0])];
					
					# if desired to overwrite input phase with GW phase, do it here
					if args.gw_phase_vcf == 1 and gw_stat >= args.gw_phase_vcf_min_confidence:
						if "-" not in gw_phase_out:
							xfields = vcf_columns[sample_column].split(":");
							new_phase = "|".join(gw_phase_out);
							if xfields[gt_index] != new_phase: phase_corrections += 1;
							xfields[gt_index] = new_phase;
							vcf_columns[sample_column] = ":".join(xfields);
						
					vcf_columns[sample_column] += ":"+"|".join(alleles_out)+":"+list_to_string(variants_out)+":"+str(block_index)+":"+"|".join(gw_phase_out)+":"+str(gw_stat);
				else:
					vcf_columns[sample_column] += ":"+"/".join(sorted(genotype))+":.:.:"+vcf_columns[sample_column].split(":")[gt_index]+":."
				
			# update other samples, that phaser was not run on, with blank values
			for i in range(9,len(vcf_columns)):
				if i != sample_column:
					vcf_columns[i] += ":.:.:.:.:.";
			
			vcf_out.write("\t".join(vcf_columns)+"\n");
					
	vcf_out.close();
	
	fun_flush_print("     corrected genome wide phase of %d variants"%(phase_corrections));	
	
def str_join(joiner,list):
	list = map(str, list);
	return(joiner.join(list));
	
def build_haplotypes(dict_variant_overlap):
	
	block_haplotypes = [];
	total_hap_pool = len(dict_variant_overlap);
	remaining_hap_pool = dict_variant_overlap;
	set_remaining_hap_pool = set(dict_variant_overlap.keys());
	while len(remaining_hap_pool) > 0:
		# this will construct many iterations of the same haplotype need to filter it out;
		# start the process with a variant pair;
		seed_var = remaining_hap_pool.keys()[0];
		seed = set([seed_var] + list(remaining_hap_pool[seed_var]));
		del remaining_hap_pool[seed_var];
		set_remaining_hap_pool.remove(seed_var);
		result = build_haplotype_v3(seed,remaining_hap_pool,set_remaining_hap_pool);
		new_hap = list(result[0]);
		# sort by location
		new_hap = sort_var_ids(new_hap);
		remaining_hap_pool = result[1];
		set_remaining_hap_pool = result[2];
		block_haplotypes.append(new_hap);
	
	return(block_haplotypes);

def sort_var_ids(ids):
	xsplit = [x.split("_") for x in ids];
	xsort = sorted(xsplit, key = lambda x: (x[0], int(x[1])))
	return(["_".join(x) for x in xsort]);

def count_hap_junctions(block):
	counted = set([]);
	reads = [];
	
	for var_index in range(0,len(block)):
		for var_allele in range(0,2):
			for other_index in range(0, len(block)):
				if other_index != var_index:
					for other_allele in range(0,2):
						if (str(var_index)+":"+str(var_allele)+":"+str(other_index)+":"+str(other_allele) not in counted) and (str(other_index)+":"+str(other_allele)+":"+str(var_index)+":"+str(var_allele) not in counted):
							reads += list(dict_variant_reads[block[var_index]]['read_set'][var_allele] & dict_variant_reads[block[other_index]]['read_set'][other_allele]);
							counted.add(str(var_index)+":"+str(var_allele)+":"+str(other_index)+":"+str(other_allele));
	return([block,len(reads)]);

def count_hap_reads(input):
	block = input[0];
	configuration = input[1];
	parent_block = None;
	block_number = None;
	
	if len(input) == 4:
		parent_block = input[2];
		block_number = input[3];
		
	global dict_variant_reads;
	
	reads = [];
	counted = set([]);
	# sum up supporting reads between all configs
	for var_index in range(0,len(block)):
		for other_index in range(0, len(block)):
			if other_index != var_index:
				if (str(var_index)+":"+str(other_index) not in counted) and (str(other_index)+":"+str(var_index) not in counted):
					# the noise test should be done here, only pairs where the signal is above noise should be counted.
					reads += list(dict_variant_reads[block[var_index]]['read_set'][int(configuration[var_index])] & dict_variant_reads[block[other_index]]['read_set'][int(configuration[other_index])]);
					counted.add(str(var_index)+":"+str(other_index));
					
	return([block, configuration, len(reads), parent_block, block_number]);

def generate_hap_network_all(input):
	block = input;
	
	global dict_variant_reads;
	
	reads = [];
	counted = set([]);
	
	out_junctions = [];
	
	for var_index in range(0,len(block)):
		for other_index in range(0, len(block)):
			if other_index != var_index:
				for allele_index in range (0,2):
					for other_allele_index in range(0,2):
						if (str(var_index)+":"+str(allele_index)+":"+str(other_index)+":"+str(other_allele_index) not in counted) and (str(other_index)+":"+str(other_allele_index)+":"+str(var_index)+":"+str(allele_index) not in counted):
							junctions = list(dict_variant_reads[block[var_index]]['read_set'][allele_index] & dict_variant_reads[block[other_index]]['read_set'][other_allele_index]);
							out_junctions.append([dict_variant_reads[block[var_index]]['id']+":"+dict_variant_reads[block[var_index]]['alleles'][allele_index],dict_variant_reads[block[other_index]]['id']+":"+dict_variant_reads[block[other_index]]['alleles'][other_allele_index], len(junctions), 0]);
							out_junctions.append([dict_variant_reads[block[var_index]]['id']+":"+dict_variant_reads[block[var_index]]['alleles'][int(not allele_index)],dict_variant_reads[block[other_index]]['id']+":"+dict_variant_reads[block[other_index]]['alleles'][int(not other_allele_index)], len(junctions), 1]);
							counted.add(str(var_index)+":"+str(allele_index)+":"+str(other_index)+":"+str(other_allele_index));
	
	return([out_junctions, block]);
	
def generate_hap_network(input):
	block = input[0];
	configuration = input[1];
	
	global dict_variant_reads;
	
	reads = [];
	counted = set([]);
	
	out_junctions = [];
	
	# sum up supporting reads between all configs
	for var_index in range(0,len(block)):
		for other_index in range(0, len(block)):
			if other_index != var_index:
				if (str(var_index)+":"+str(other_index) not in counted) and (str(other_index)+":"+str(var_index) not in counted):
					
					## SHOULD FIRST CHECK TO MAKE SURE THIS ISN'T A READ PAIR THAT FAILED THE TEST
					# actually I don't think this matters, it will always choose the most supported phase
					
					junctions = list(dict_variant_reads[block[var_index]]['read_set'][int(configuration[var_index])] & dict_variant_reads[block[other_index]]['read_set'][int(configuration[other_index])]);
					out_junctions.append([dict_variant_reads[block[var_index]]['rsid']+":"+dict_variant_reads[block[var_index]]['alleles'][int(configuration[var_index])],dict_variant_reads[block[other_index]]['rsid']+":"+dict_variant_reads[block[other_index]]['alleles'][int(configuration[other_index])], len(junctions), 0]);
					out_junctions.append([dict_variant_reads[block[var_index]]['rsid']+":"+dict_variant_reads[block[var_index]]['alleles'][int(not int(configuration[var_index]))],dict_variant_reads[block[other_index]]['rsid']+":"+dict_variant_reads[block[other_index]]['alleles'][int(not int(configuration[other_index]))], len(junctions), 1]);
					counted.add(str(var_index)+":"+str(other_index));
					
	return([out_junctions, block, configuration]);
	
def get_allele_phase(allele, var_dict):
	
	try:
		return(var_dict['phase'].index(allele));
	except:
		return(float('nan'));

def build_haplotype_v3(set_haplotype, dict_all_associations, set_all_associations):
	global args;
	
	overlapping = set_haplotype & set_all_associations;
	
	while len(overlapping) > 0:
		for variant in overlapping:
			set_haplotype = set_haplotype | dict_all_associations[variant];
			del dict_all_associations[variant];
			set_all_associations.remove(variant);
		
		overlapping = set_haplotype & set_all_associations;
		
	return([set_haplotype, dict_all_associations, set_all_associations])

def get_var_pos(var_fields):
	return(int(var_fields.split(":")[0]));

def list_to_string(xlist,sep=","):
	string_out = "";
	for item in xlist:
		string_out += str(item) + sep;
	
	if len(sep) > 0:
		string_out = string_out[:-len(sep)];
	
	return(string_out);

def intersect_reads_vars_file(file):
	stream_in = open(file, "r");
	output = intersect_reads_vars(stream_in);
	stream_in.close();
	return(output);

def intersect_reads_vars(reads):
	global args;
	global dict_vcf_interval;
	global dict_variant_reads;
	global as_cutoff;
	global use_as_cutoff;
	global ab_cutoff;
	global isize_cutoff;
	
	out_intersect = [];
	out_orphans = "";
	
	#these are used to estimate sequencing noise level
	base_matches = {};
	base_mismatches = {};
	variant_overlaps = [];
	
	chr = "";
	map_start = 0;
	
	for read in reads:
		read_columns = read.replace("\n","").split("\t");
		if len(read_columns) > 1:
			#ERR188213.16580131	419	1	77078	0	20M242102N56M	=	319229	242227	GGTGGGAGGATCGCTTGAACTCAGGAGTTTGAGACCA
			read_id = read_columns[0];
			chr = read_columns[2];
			
			map_start = int(read_columns[3]);
			cigar = read_columns[5]
			template_length = abs(int(read_columns[8]));
			
			# get the AS
			alignment_score = 0;
			
			for i in range(11,len(read_columns)):
				if "AS:" in read_columns[i]:
					alignment_score = int(read_columns[i].split(":")[2]);
			
			variant_overlap = 0;
			
			# make sure the read passes isize and AS cutoffs
			if (isize_cutoff == 0 or template_length <= isize_cutoff) and (use_as_cutoff == False or alignment_score >= as_cutoff):
				#process the cigar string and use to build a pseudoread
				number_build = "";
				read_pos = 0;
				
				#first filter read by BASEQ
				read_seq = "";
				for base, quality in zip(read_columns[9], read_columns[10]):
					baseq = ord(quality) - 33;
					if baseq >= args.baseq:
						read_seq += base;
					else:
						read_seq += "N";
					
				aligned_bases = 0;
				pseudo_read = "";
				insertions = {};
				for c in cigar:
					ascii = ord(c);
					if ascii >= 48 and ascii <= 57:
						#number;
						number_build += c;
					else:
						seq_len = int(number_build);
						if c == "M":
							#match
							pseudo_read += read_columns[9][read_pos:read_pos+seq_len];
							read_pos += seq_len;
							aligned_bases += seq_len;
						elif c == "N":
							#gap
							pseudo_read += "N" * seq_len;
						elif c == "D":
							#deletion
							pseudo_read += "D" * seq_len;
						elif c == "I":
							#insertion
							insertions[read_pos] = read_columns[9][read_pos:read_pos+seq_len];
							read_pos += seq_len;
						elif c == "S":
							#softclip, reads not used in alignment
							read_pos += seq_len;
						
				
						number_build = "";
				# make sure the number of aligned bases meets minimum requirements
				if ab_cutoff == float('nan') or aligned_bases >= ab_cutoff:				
					#now get all the variants that the read overlaps
					overlapped_variants = [];
					if chr in dict_vcf_interval:							
						variants = dict_vcf_interval[chr][map_start-1:map_start+len(pseudo_read)-1];
							
						#for each variant check the base of the read and assign to a haplotype
				
						for interval in variants:
							var_id = interval.data;							
							var_dict = dict_variant_reads[var_id];
							var_seq = var_dict['alleles'];
							read_start = (interval.begin+1)-map_start;
							read_end = (interval.end+1)-map_start;
						
							try:
								read_seq = pseudo_read[read_start:read_end]
							except:
								read_seq = "";
							
							for read_pos in range(read_start,read_end+1):
								if read_pos in insertions:
									read_seq += insertions[read_pos];
							
							#remove deletion placeholder
							read_seq = read_seq.replace("D","");
							if read_seq != "N":
								#update the dictionary with reads that match this allele
								if read_seq in var_seq:
									for i in range(0,len(var_seq)):
										if var_seq[i] == read_seq:
											out_intersect.append([var_id,i,read_id]);
											if var_id not in base_matches: base_matches[var_id] = 0;
											base_matches[var_id] += 1;
											variant_overlap += 1;
								else:
									#read doesn't map to either allele
									#record this to the other category
									out_intersect.append([var_id,-1,read_id]);
									if var_id not in base_mismatches: base_mismatches[var_id] = 0;
									base_mismatches[var_id] += 1;
									if args.output_orphans == 1:
										out_orphans += read+"\n";
					
						if variant_overlap > 1:
							variant_overlaps.append(variant_overlap);
					
	
	print_debug("     completed %s - %d"%(chr,map_start));
	
	return([out_intersect,out_orphans,base_matches,base_mismatches,variant_overlaps]);

def print_warning(text):
	if args.show_warning == 1:
		fun_flush_print(text);

def dict_from_info(info_field):
	out_dict = {};
	fields = info_field.split(";");
	for field in fields:
		sub_field = field.split("=");
		if len(sub_field) == 2:
			out_dict[sub_field[0]] = sub_field[1];
	
	return(out_dict);

def fun_flush_print(text):
	print(text);
	sys.stdout.flush();

def fatal_error(text):
	fun_flush_print("     FATAL ERROR: "+text);
	quit();

def print_debug(text):
	if args.debug == 1:
		fun_flush_print(text);

def pool_split(threads, data):
	global args;
	data_length = len(data);
	pool_input = [];
	
	# calculate pool size if all data is divided by number of threads
	optimal_pool_size = data_length/threads;
	
	# unfortunately due to OS limitations the maximum output but a given thread is limited
	# so the pool size can't be too enormous
	# see : http://bugs.python.org/issue8426
	# so limit it to at max some value (set at 100,000 by default)
	# this is probably conservative but I haven't checked out what the best number is yet
	
	pool_size = min([args.max_items_per_thread, optimal_pool_size]);
	pool_inputs = data_length / pool_size;
	
	for i in range(0,pool_inputs):
		#last pool gets the remaining reads
		if i == (pool_inputs-1):
			pool_input.append(data[(i*pool_size):]);
		else:
			pool_input.append(data[(i*pool_size):((i+1)*pool_size)]);
	
	return(pool_input);

def pool_setup(pool_input):
	global args;
	
	threads = min([len(pool_input),args.threads]);
	
	return (multiprocessing.Pool(processes=threads));

def parallelize(function, pool_input):
	global args;
	
	threads = min([len(pool_input),args.threads]);
	if args.threads > 1:
		pool = multiprocessing.Pool(processes=threads);
		pool_output = pool.map(function, pool_input);
		pool.close() # no more tasks
		pool.join()  # wrap up current tasks
	else:
		pool_output = [];
		for input in pool_input:
			pool_output.append(function(input));
	
	return(pool_output);

def generate_intervals(input):
	chr = input[0];
	out_interval = IntervalTree();
	
	for variant in input[1:]:
		id = variant[0];
		pos = int(variant[1]);
		length = int(variant[2]);
		out_interval[pos-1:pos+(length-1)] = id;
	
	print_debug("          completed interval for chr: %s"%(chr));
	
	return([chr, out_interval]);
	
def process_vcf(lines):
	global dict_blacklist_interval;
	global args;
	global sample_column;
	global vcf_annot;
	
	data_out = [];
	black_listed = 0;
	total_hets = 0;
	indels_excluded = 0;
	
	for line in lines:
		if line[0:1] != "#":
			##CHROM  POS     ID      REF     ALT     QUAL    FILTER  INFO    FORMAT  NA06986
			vcf_columns = line.replace("\n","").split("\t");
			if len(vcf_columns) > sample_column:
				id = vcf_columns[2];
				chrom = vcf_columns[0];
				pos = int(vcf_columns[1]);
			
				if (args.chr == "" or chrom == args.chr) and (vcf_columns[6] == "PASS" or args.pass_only == 0):
					if "GT" in vcf_columns[8]:
						sample_columns = vcf_columns[8].split(":");
						gt_index = sample_columns.index("GT");
						genotype = list(vcf_columns[sample_column].split(":")[gt_index]);
						
						is_phased = 0;
						if "|" in genotype:
							is_phased = 1;
							genotype.remove("|");
						if "/" in genotype: genotype.remove("/");
						
						info_fields = annotation_to_dict(vcf_columns[7])
						maf = None;
						if args.gw_af_field in info_fields:
							# make sure to get the right index if multi-allelic site
							afs = map(float, info_fields[args.gw_af_field].split(","));
							use_afs = [];
							for allele in list(genotype):
								if allele != "." and int(allele) != 0:
									use_afs.append(int(allele) - 1);
							# if there are multiple alternative alleles use the lowest MAF
							if len(use_afs) > 0:
								maf = min([min([afs[x],1-afs[x]]) for x in use_afs]);
						
						if args.only_phased == 0 or is_phased == 1:
							if len(set(genotype)) > 1:
								#this is a het site;
								total_hets += 1;
								
								#check to see if SNP falls into a blacklisted region
								if chrom in dict_blacklist_interval:
									blacklist = dict_blacklist_interval[chrom][pos-1:pos];
								else:
									blacklist = [];
		
								if len(blacklist) == 0:
									if chrom not in dict_vcf_interval: dict_vcf_interval[chrom] = IntervalTree();

									# get only the alleles this individual has
									alt_alleles = vcf_columns[4].split(",");
									all_alleles = [vcf_columns[3]] + alt_alleles;
									ind_alleles = [];

									for i in range(0,len(all_alleles)):
										if str(i) in genotype:
											ind_alleles.append(all_alleles[i]);
					
									if len(ind_alleles) == 2:					
										phase = [];
										if is_phased == 1:
											for index in genotype:
												phase.append(all_alleles[int(index)]);
										else:
											phase = ["-","-"];
					
										#get length of reference sequence
										ref_length = len(vcf_columns[3]);
										all_lengths = [len(x) for x in ind_alleles];
										if max(all_lengths) == 1 or args.include_indels == 1:
											#generate a unique id
											unique_id = chrom+"_"+str(pos)+"_"+("_".join(all_alleles));
						
											data_out.append({"id":unique_id, "rsid":id,"ref":vcf_columns[3],"chr":chrom,"pos":pos,"length":ref_length,"alleles":ind_alleles,"phase":phase, "gw_phase":phase, "maf":maf, "other_reads":[], "reads":[[] for i in range(len(ind_alleles))]});
										else:
											indels_excluded += 1;
																		
								else:
									black_listed += 1;
									print_warning("          warning: blacklisted variant %s - %s"%(id,list(blacklist)[0].data));
						else:
							print_warning("          warning: unphased variant %s"%(chrom+":"+str(pos)+"-"+id));
	
	return([total_hets, black_listed, indels_excluded, data_out]);

def calculate_as_cutoff(reads):
	global args;
	
	if args.as_q_cutoff > 0:
		alignment_scores = [];
	
		for read in reads:
			read_columns = read.replace("\n","").split("\t");
			if len(read_columns) > 1:
				as_found = 0;
				for i in range(11,len(read_columns)):
					if "AS:" in read_columns[i]:
						value = read_columns[i].split(":")[2];
						alignment_scores.append(int(value));
						as_found = 1;
				if as_found == 0:
					# there is no alignment score in the file
					fun_flush_print("          no alginment score found in BAM file, can't use AS cutoff.");
					return([False,0]);
				
		percentile_cutoff = numpy.percentile(alignment_scores,args.as_q_cutoff*100);
	
		#stream_out.close();
		return([True,percentile_cutoff]);
	else:
		return([False,0]);

def calculate_ab_cutoff(reads):
	global args;
	
	if args.ab_q_cutoff > 0:
		aligned_bases = [];
		for read in reads:
			read_columns = read.replace("\n","").split("\t");
			if len(read_columns) > 1:
				read_ab = 0;
				#ERR188213.16580131	419	1	77078	0	20M242102N56M	=	319229	242227	GGTGGGAGGATCGCTTGAACTCAGGAGTTTGAGACCA
				cigar = read_columns[5];
				number_build = "";
				for c in cigar:
					ascii = ord(c);
					if ascii >= 48 and ascii <= 57:
						#number;
						number_build += c;
					else:
						seq_len = int(number_build);
						if c == "M":
							#match
							read_ab += seq_len;
						
						number_build = "";
				
				aligned_bases.append(read_ab);
		
		percentile_cutoff = numpy.percentile(aligned_bases,args.ab_q_cutoff*100);
		return(percentile_cutoff);
	else:
		return(0);
		
def annotation_to_dict(text,sep=";"):
	dict_out = {};
	vars = text.split(sep);
	for var in vars:
		if "=" in var:
			key = var.split("=")[0];
			values = var.split("=")[1];
		dict_out[key] = values;
	return(dict_out);
								
if __name__ == "__main__":
	main();