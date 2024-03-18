import pysam;
import argparse;
import vcf;
import gzip;
import copy;
import multiprocessing;
import subprocess;
import sys;

def main():
	global vcf_af;
	global args;

	parser = argparse.ArgumentParser()
	# required
	parser.add_argument("--geno_vcf", help="VCF containing phased genotype.")
	parser.add_argument("--sample", help="Name of sample to use in VCF file.")
	parser.add_argument("--af_vcf", help="VCF to retrieve allele frequencies from. Must be indexed with tabix. If left blank will attempt to retrieve from genotype vcf.")
	parser.add_argument("--af_field", default='AF', help="Allele frequency field in af_vcf to use ('AF' by default).")
	parser.add_argument("--cadd_file", help="The path to the CADD 'whole_genome_SNVs.tsv.gz' file, which contains a CADD score and annotation for every SNV in the genome.")
	parser.add_argument("--o", help="Output file")
	parser.add_argument("--threads", type=int, default=1, help="Number of threads to use.")
	args = parser.parse_args()

	# initialize AF vcf
	if args.af_vcf != None:
		vcf_af = vcf.Reader(filename=args.af_vcf);

	if args.o == None:
		print("Error: please specify an output directory.");
		sys.exit(1)

	# first assign each variant to a gene and record CADD info
	print("1. Reading VCF...");
	if "gz" in args.geno_vcf:
		stream_vcf = gzip.open(args.geno_vcf,"r");
	else:
		stream_vcf = open(args.geno_vcf,"r");

	global dict_gw_variant_info;
	global dict_pg_variant_info;
	global dict_gw_gene_variants;
	global dict_pg_gene_variants;
	global dict_uniqueid_rsid;

	dict_gw_variant_info = {};
	dict_pg_variant_info = {};
	dict_gw_gene_variants = {};
	dict_pg_gene_variants = {};
	dict_uniqueid_rsid = {};

	cadd_retrieve_list_gw = [];
	cadd_retrieve_list_pg = [];

	sample_column = 0;
	for line in stream_vcf:
		columns = line.replace("\n","").split("\t");

		if line[0:4] == "#CHR":
			if args.sample in columns:
				sample_column = columns.index(args.sample);
			else:
				print("Error sample not found in VCF.")
				sys.exit(1)
		elif line[0:1] != "#":
			# lookup each heterozygous variant
			chrom = columns[0];
			pos = columns[1];
			rsid = columns[2];
			ref = columns[3];
			alt = columns[4];
			unique_id = "_".join([chrom,pos,ref,alt]);
			info_fields = annotation_to_dict(columns[7])

			dict_uniqueid_rsid[unique_id] = rsid;

			if len(columns[8].split(":")) == len(columns[sample_column].split(":")):
				# use both the genome wide phase
				gt_index = columns[8].split(":").index("GT");
				genotype = list(columns[sample_column].split(":")[gt_index]);

				# require a full genotype
				# done care about homo refs
				if "." not in genotype and genotype.count("0") != 2:
					gt_alleles = copy.deepcopy(genotype);
					if "/" in gt_alleles: gt_alleles.remove("/");
					if "|" in gt_alleles: gt_alleles.remove("|");

					if "|" in genotype or len(set(gt_alleles)) == 1:
						# individual is either homozygote or phased
						# now retrieve the CADD annotation
						cadd_retrieve_list_gw.append([unique_id,info_fields,gt_alleles,0]);

				# and the phaser phase
				if "PG" in columns[8].split(":"):
					gt_index = columns[8].split(":").index("PG");
					genotype = list(columns[sample_column].split(":")[gt_index]);

					if "." not in genotype and genotype.count("0") != 2 and "/" not in genotype:
						pi_index = columns[8].split(":").index("PI");
						block_index = float(columns[sample_column].split(":")[pi_index]);

						gt_alleles = copy.deepcopy(genotype);
						if "|" in gt_alleles: gt_alleles.remove("|");

						if "|" in genotype or len(set(gt_alleles)) == 1:
							# individual is either homozygote or phased
							# now retrieve the CADD annotation
							cadd_retrieve_list_pg.append([unique_id,info_fields,gt_alleles,block_index]);

			else:
				print(("Column info error %s"%(unique_id)));

	stream_vcf.close();

	print("2. Retrieving CADD info for all phased variants...");
	# first retrieve GW

	pool = multiprocessing.Pool(processes=args.threads);
	pool_output = pool.map(get_variant_cadd, cadd_retrieve_list_gw);
	pool.close() # no more tasks
	pool.join()  # wrap up current tasks

	for result in pool_output:
		unique_id = result[0];
		gt_alleles = result[1];
		cadd_info = result[2];
		gene_list = result[3];
		phaser_bi = result[4];

		dict_gw_variant_info[unique_id] = [gt_alleles] + [cadd_info] + [gene_list] + [phaser_bi];

		for gene in gene_list:
			if gene not in dict_gw_gene_variants: dict_gw_gene_variants[gene] = [];
			dict_gw_gene_variants[gene].append(unique_id);


	set_retrieved_vars = set(dict_gw_variant_info.keys());

	cadd_retrieve_list_pg_flt = [];

	# don't retrieve the info twice
	# check if we have already retrieved or not
	for variant in cadd_retrieve_list_pg:
		unique_id = variant[0];
		if unique_id in set_retrieved_vars:
			dict_pg_variant_info[unique_id] = dict_gw_variant_info[unique_id];
			cadd_info = dict_gw_variant_info[unique_id][1];
			gene_list = dict_gw_variant_info[unique_id][2];

			for gene in gene_list:
				if gene not in dict_pg_gene_variants: dict_pg_gene_variants[gene] = [];
				dict_pg_gene_variants[gene].append(unique_id);
		else:
			cadd_retrieve_list_pg_flt.append(variant);

	if len(cadd_retrieve_list_pg_flt) > 0:
		# now retrieve remaining PG variants
		pool = multiprocessing.Pool(processes=args.threads);
		pool_output = pool.map(get_variant_cadd, cadd_retrieve_list_pg_flt);
		pool.close() # no more tasks
		pool.join()  # wrap up current tasks

		for result in pool_output:
			unique_id = result[0];
			gt_alleles = result[1];
			cadd_info = result[2];
			gene_list = result[3];
			phaser_bi = result[4];

			dict_pg_variant_info[unique_id] = [gt_alleles] + [cadd_info] + [gene_list] + [phaser_bi];
			for gene in gene_list:
				if gene not in dict_pg_gene_variants: dict_pg_gene_variants[gene] = [];
				dict_pg_gene_variants[gene].append(unique_id);

	print("3. Retrieving variant allele frequencies...");
	global dict_allele_af;
	dict_allele_af = {};

	if args.af_vcf != None:
		# first build the list of allele freqs to get
		af_retrieve_list = set([]);

		for variant in list(dict_gw_variant_info.keys()):
			for allele in dict_gw_variant_info[variant][1]:
				chr = variant.split("_")[0];
				pos = variant.split("_")[1];
				af_retrieve_list.add(chr+"_"+pos+"_"+dict_gw_variant_info[variant][1][allele][7]);

		for variant in list(dict_pg_variant_info.keys()):
			for allele in dict_pg_variant_info[variant][1]:
				chr = variant.split("_")[0];
				pos = variant.split("_")[1];
				af_retrieve_list.add(chr+"_"+pos+"_"+dict_pg_variant_info[variant][1][allele][7]);

		if len(af_retrieve_list) > 0:
			pool = multiprocessing.Pool(processes=args.threads);
			pool_output = pool.map(get_variant_af, af_retrieve_list);
			pool.close() # no more tasks
			pool.join()  # wrap up current tasks

			for result in pool_output:
				dict_allele_af[result[0]] = result[1];

	# now
	print("4. Identifying cases of compound heterozygosity...");
	all_genes = set(dict_gw_gene_variants.keys()) | set(dict_pg_gene_variants.keys());

	# not meaningful to find interactions in unannotated variants
	all_genes.remove("NA");

	pool = multiprocessing.Pool(processes=args.threads);
	pool_output = pool.map(get_gene_interactions, all_genes);
	pool.close() # no more tasks
	pool.join()  # wrap up current tasks

	stream_out = open(args.o, "w");
	stream_out.write("\t".join(["ensg","name","variant_a","rsid_a","allele_a","af_a","cadd_phred_a","cadd_effect_a","variant_b","rsid_b","allele_b","af_b","cadd_phred_b","cadd_effect_b","configuration","read_backed"])+"\n");

	for gene_result in pool_output:
		for interaction in gene_result:
			stream_out.write("\t".join(map(str,interaction))+"\n");
	stream_out.close();

def get_interactions(variant_a,variant_b):
	out_interactions = [];

	# first check to see if they are in the same phase block
	# note if this is coming from genome wide phasing it is always set to 0
	# so they will be on the same block
	# otherwise it is set to the phaser block index

	if variant_a[3] == variant_b[3]:
		for index_a in range(0,len(variant_a[0])):
			for index_b in range(0,len(variant_a[0])):
				if index_a == index_b:
					out_interactions.append([int(variant_a[0][index_a]),int(variant_b[0][index_b]),"cis"]);
				else:
					out_interactions.append([int(variant_a[0][index_a]),int(variant_b[0][index_b]),"trans"]);

	# don't care about interactions involving reference
	keep_interactions = [];
	for interaction in out_interactions:
		if interaction[0] != 0 and interaction[1] != 0:
			keep_interactions.append(interaction);

	return(keep_interactions);

def get_variant_af(info):
	global vcf_af;
	global args;

	chr = info.split("_")[0];
	pos = int(info.split("_")[1])
	allele = info.split("_")[2];

	key = chr+"_"+str(pos)+"_"+allele;
	# retrieve AF
	variant_afs = [];
	records = vcf_af.fetch(chr,pos-1,pos);

	for record in records:
		if record.POS == pos:
			alts = record.ALT;
			variant_afs = record.INFO[args.af_field];
			# in some VCFs this is returned as a list, in others not...
			if isinstance(variant_afs, list) == False:
				variant_afs = [afs];
			break;

	if len(variant_afs) == 0:
		return([key,0]);
	else:
		if allele in alts:
			allele_index = alts.index(allele);
			allele_freq = variant_afs[allele_index];
			# don't want maf, should be the AF for this particular allele.
			#maf = min([1-allele_freq,allele_freq]);
			return([key,allele_freq]);
		else:
			return([key,0]);

def get_variant_cadd(input):
	global args;

	unique_id,info_fields,gt_alleles,phaser_bi = input;

	variant = unique_id.split("_");
	chr = variant[0];
	pos = int(variant[1]);
	alleles = [variant[2],variant[3]];
	alt_alleles = alleles[1].split(",");

	output = {};
	gene_list = [];

	# this turns out to be slower than the solution below
	#records = subprocess.check_output("set -euo pipefail && "+"tabix "+args.cadd_file+" "+chr+":"+str(pos-1)+"-"+str(pos), shell=True, executable='/bin/bash').split("\n");

	# initialize CADD retrieval
	# have to do this for every variant because of problem with pysam
	tabix_cadd = pysam.Tabixfile(args.cadd_file,"r");
	records = tabix_cadd.fetch(chr,pos-1,pos);

	#Chrom  Pos     Ref     Anc     Alt     Type    Length  isTv    isDerived       AnnoType        Consequence     ConsScore       ConsDetail      GC      CpG     mapAbility20bp  mapAbility35bp  scoreSegDup     priPhCons       mamPhCons       verPhCons       priPhyloP       mamPhyloP       verPhyloP       GerpN   GerpS   GerpRS  GerpRSpval      bStatistic      mutIndex        dnaHelT dnaMGW  dnaProT dnaRoll mirSVR-Score    mirSVR-E        mirSVR-Aln      targetScan      fitCons cHmmTssA        cHmmTssAFlnk    cHmmTxFlnk      cHmmTx  cHmmTxWk        cHmmEnhG        cHmmEnh cHmmZnfRpts     cHmmHet cHmmTssBiv      cHmmBivFlnk     cHmmEnhBiv      cHmmReprPC      cHmmReprPCWk    cHmmQuies       EncExp  EncH3K27Ac      EncH3K4Me1      EncH3K4Me3      EncNucleo       EncOCC  EncOCCombPVal   EncOCDNasePVal  EncOCFairePVal  EncOCpolIIPVal  EncOCctcfPVal   EncOCmycPVal    EncOCDNaseSig   EncOCFaireSig   EncOCpolIISig   EncOCctcfSig    EncOCmycSig     Segway  tOverlapMotifs  motifDist       motifECount     motifEName      motifEHIPos     motifEScoreChng TFBS    TFBSPeaks       TFBSPeaksMax    isKnownVariant  ESP_AF  ESP_AFR ESP_EUR TG_AF   TG_ASN  TG_AMR  TG_AFR  TG_EUR  minDistTSS      minDistTSE      GeneID  FeatureID       CCDS    GeneName        cDNApos relcDNApos      CDSpos  relCDSpos       protPos relProtPos      Domain  Dst2Splice      Dst2SplType     Exon    Intron  oAA     nAA     Grantham        PolyPhenCat     PolyPhenVal     SIFTcat SIFTval RawScore        PHRED
	#1       69139   C       NA      T       SNV     0       FALSE   NA      CodingTranscript        STOP_GAINED     8       stop_gained     0.37    0.04    0.333333        0.333333        0.99    0.020   0.961   0.109   0.393   1.228   1.785   2.31    2.31    928     5.41355e-97     994     -4      1.29    0.62    -2.97   1.20    NA      NA      NA      NA      0.497099        0.000   0.000   0.000   0.000   0.000   0.000   0.000   0.000   0.000   0.000   0.000   0.000   0.000   0.000   1.000   0.41    25.60   32.04   21.00   NA      NA      NA      NA      NA      NA      NA      NA      NA      NA      NA      NA      NA      D       NA      NA      NA      NA      NA      NA      NA      NA      NA      FALSE   NA      NA      NA      NA      NA      NA      NA      NA      48      870     ENSG00000186092 ENST00000335137 CCDS30547.1     OR4F5   49      0.05    49
	for record in records:
		# NOTE THIS WILL ONLY GET THE ANNOTATION FOR THE ALTERNATIVE ALLELE
		if record != "":
			vfields = record.rstrip().split('\t');
			if int(vfields[1]) == pos:
				if vfields[4] in alt_alleles:
					a_index = alt_alleles.index(vfields[4]) + 1;
					phred = vfields[len(vfields)-1];
					annotation = vfields[10];
					gene_ensg = vfields[92];
					gene_name = vfields[95];
					var_alt = vfields[4];

					gene_list.append(gene_ensg);

					allele_freq = None;
					#if possible get maf from input VCF
					if args.af_vcf == None:
						if args.af_field in info_fields:
							# make sure to get the right index if multi-allelic site
							afs = list(map(float, info_fields[args.af_field].split(",")));
							allele_freq = afs[a_index-1];
							#maf = min([1-allele_freq,allele_freq]);

					output[gene_ensg+":"+str(a_index)] = [phred,annotation,gene_ensg,gene_name,chr,pos,allele_freq,var_alt];
	return([unique_id,gt_alleles,output,gene_list,phaser_bi]);

def annotation_to_dict(text,sep=";"):
	dict_out = {};
	vars = text.split(sep);
	for var in vars:
		if "=" in var:
			key = var.split("=")[0];
			values = var.split("=")[1];
			dict_out[key] = values;
	return(dict_out);

def get_gene_interactions(xgene):
	global dict_gw_variant_info;
	global dict_pg_variant_info;
	global dict_gw_gene_variants;
	global dict_pg_gene_variants;
	global dict_uniqueid_rsid;
	global dict_allele_af;

	out_interactions = [];
	outputted_configs = set([]);

	# first get compound hets from GW phasing
	# also output any read backed phasing which disagrees
	if xgene in dict_gw_gene_variants:
		for variant_a in dict_gw_gene_variants[xgene]:
			for variant_b in dict_gw_gene_variants[xgene]:
				if variant_a != variant_b:
					# 1 determine type of interaction
					gw_interactions = get_interactions(dict_gw_variant_info[variant_a],dict_gw_variant_info[variant_b]);

					# 2 determine if read backed supports or disagrees with GW
					rb_interactions = [];
					if xgene in dict_pg_gene_variants:
						if variant_a in dict_pg_gene_variants[xgene] and variant_b in dict_pg_gene_variants[xgene]:
							rb_interactions = get_interactions(dict_pg_variant_info[variant_a],dict_pg_variant_info[variant_b]);

					read_backed = "0";
					if len(gw_interactions) == len(rb_interactions) and gw_interactions == rb_interactions:
						# read backed supported and concordant
						read_backed = "1";
					if len(gw_interactions) == len(rb_interactions) and gw_interactions != rb_interactions:
						# read backed data available, but conflicting
						read_backed = "-1";
					elif len(rb_interactions) == 0:
						# not read backed
						read_backed = "0";

					out_interactions += build_interaction_result(xgene, variant_a, dict_gw_variant_info[variant_a], variant_b, dict_gw_variant_info[variant_b],gw_interactions,read_backed);

					if read_backed == "-1":
						out_interactions += build_interaction_result(xgene, variant_a, dict_gw_variant_info[variant_a], variant_b, dict_gw_variant_info[variant_b],rb_interactions,"1");

					outputted_configs.add(variant_a+"_"+variant_b);

	if xgene in dict_pg_gene_variants:
		# now get interactions from phASER phasing, that were not previously outputted

		for variant_a in dict_pg_gene_variants[xgene]:
			for variant_b in dict_pg_gene_variants[xgene]:
				if variant_a != variant_b:
					# only output interactions not previously outputted in GW
					if variant_a+"_"+variant_b not in outputted_configs:
						# 1 get interactions
						pg_interactions = get_interactions(dict_pg_variant_info[variant_a],dict_pg_variant_info[variant_b]);
						# build result
						out_interactions += build_interaction_result(xgene, variant_a, dict_pg_variant_info[variant_a], variant_b, dict_pg_variant_info[variant_b],pg_interactions,"1");

						outputted_configs.add(variant_a+"_"+variant_b);

	return(out_interactions);

def build_interaction_result(xgene, variant_a, dict_variant_info_a, variant_b, dict_variant_info_b, interactions, read_backed):
	global dict_allele_af;

	output = [];

	for interaction in interactions:
		allele_a = interaction[0];
		allele_b = interaction[1];
		configuration = interaction[2];

		# need to ensure that there is an annotation for both alleles
		# otherwise this won't work
		# in some cases there will not be an annotation for both alleles
		# for example
		#       * allele_a
		#	----------------------------------- Gene X
		#	----------------- Gene Y
		#                             ^ allele_b
		# actually in this case for gene Y: allele_b won't have been in the same gene variant dictionary as allele_a so it shouldn't
		# ever try to get the annotations, however gene X  will have both

		if xgene+":"+str(allele_a) in dict_variant_info_a[1] and xgene+":"+str(allele_b) in dict_variant_info_b[1]:
			#[phred,annotation,gene_ensg,gene_name,chr,pos,af];
			name = dict_variant_info_a[1][xgene+":"+str(allele_a)][3];
			chr_a = dict_variant_info_a[1][xgene+":"+str(allele_a)][4];
			pos_a = dict_variant_info_a[1][xgene+":"+str(allele_a)][5];
			cadd_phred_a = dict_variant_info_a[1][xgene+":"+str(allele_a)][0];
			cadd_effect_a = dict_variant_info_a[1][xgene+":"+str(allele_a)][1];

			chr_b = dict_variant_info_b[1][xgene+":"+str(allele_b)][4];
			pos_b = dict_variant_info_b[1][xgene+":"+str(allele_b)][5];
			cadd_phred_b = dict_variant_info_b[1][xgene+":"+str(allele_b)][0];
			cadd_effect_b = dict_variant_info_b[1][xgene+":"+str(allele_b)][1];

			# get mafs
			af_a = ".";
			if args.af_vcf != None:
				allele_a = dict_variant_info_a[1][xgene+":"+str(allele_a)][7];
				af_a = dict_allele_af[chr_a+"_"+str(pos_a)+"_"+allele_a];
			elif dict_variant_info_a[1][xgene+":"+str(allele_a)][6] != None:
				af_a = dict_variant_info_a[1][xgene+":"+str(allele_a)][6];

			af_b = ".";
			if args.af_vcf != None:
				allele_b = dict_variant_info_b[1][xgene+":"+str(allele_b)][7];
				af_b = dict_allele_af[chr_b+"_"+str(pos_b)+"_"+allele_b];
			elif dict_variant_info_b[1][xgene+":"+str(allele_b)][6] != None:
				af_b = dict_variant_info_b[1][xgene+":"+str(allele_b)][6];

			output += [[xgene,name,variant_a,dict_uniqueid_rsid[variant_a],allele_a,af_a,cadd_phred_a,cadd_effect_a,variant_b,dict_uniqueid_rsid[variant_b],allele_b,af_b,cadd_phred_b,cadd_effect_b,configuration,read_backed]];

	return(output);

if __name__ == "__main__":
	main();
