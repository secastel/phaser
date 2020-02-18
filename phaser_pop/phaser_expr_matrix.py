import subprocess;
import os;
import pandas;
import multiprocessing;
import subprocess;
import argparse;

## THIS NEEDS TO BE UPDATED TO USE PROPER TEMPORARY FILE SYSTEMS

def main():
	parser = argparse.ArgumentParser()
	parser.add_argument("--gene_ae_dir", required=True, help="Directory containing phASER gene AE files.")
	parser.add_argument("--features", required=True, help="BED feature file used to generate phASER gene AE files.")
	parser.add_argument("--t", required=False, default=1, type=int, help="Number of threads to use.")
	parser.add_argument("--o", required=True, help="Output file prefix.")

	args = parser.parse_args()

	result_dir = args.gene_ae_dir;
	
	# make needed temporary directories
	if not os.path.isdir("tmp"): subprocess.call("mkdir tmp", shell=True);
	if not os.path.isdir("tmp/all"): subprocess.call("mkdir tmp/all", shell=True);
	if not os.path.isdir("tmp/gw_phased"): subprocess.call("mkdir tmp/gw_phased", shell=True);
	
	version = "0.1.0";
	print("");
	print("##################################################")
	print("        Welcome to phASER-POP v%s"%(version));
	print("  Author: Stephane Castel (scastel@nygenome.org)")
	print("##################################################");
	print("");

	print("#1 Loading gene feature file...");
	df_features = pandas.read_csv(args.features, sep="\t", index_col=False,header=None, comment="#");
	global gene_list;
	gene_list = df_features[3].tolist()

	print("#2 Loading gene ae files...");
	xinput = [];

	i = 1;

	for filename in os.listdir(result_dir):
		if ".txt" in filename:
			xinput.append([i, os.path.join(result_dir, filename)]);
			i += 1;
	
	if len(xinput) == 0:
		print("FATAL ERROR - no files read for input...");
		quit();

	xresults_raw = parallelize(read_result, args.t, xinput);

	xresults = [];
	for xresult in xresults_raw:
		for xi in xresult:
			xresults.append(xi);

	print("#3 Saving sample matrix (all)...")
	global df_matrix;
	df_matrix = generate_basic_matrix(xinput[0][1]);

	for xsample in xresults:
		xdf = pandas.read_csv("tmp/all/"+xsample+".txt", sep="\t", index_col=False);
		df_matrix[xsample] = xdf[xsample];
	df_matrix.to_csv(args.o+".bed", sep="\t", index=False);
	subprocess.call("bgzip -f "+args.o+".bed; tabix -p bed -f "+args.o+".bed.gz", shell=True);
	subprocess.call("rm tmp/all/*.txt", shell=True);

	print("#4 Saving sample matrix (gw_phased)...")
	df_matrix = generate_basic_matrix(xinput[0][1]);
	for xsample in xresults:
		xdf = pandas.read_csv("tmp/gw_phased/"+xsample+".txt", sep="\t", index_col=False);
		df_matrix[xsample] = xdf[xsample];
	df_matrix.to_csv(args.o+".gw_phased.bed", sep="\t", index=False);
	subprocess.call("bgzip -f "+args.o+".gw_phased.bed; tabix -p bed -f "+args.o+".gw_phased.bed.gz", shell=True);
	subprocess.call("rm tmp/gw_phased/*.txt", shell=True);

def generate_basic_matrix(path):
	xdf = pandas.read_csv(path, sep="\t", index_col = False);
	one_bam = list(set(xdf['bam'].tolist()))[0];
	xdf = xdf[(xdf.bam == one_bam)];
	df_matrix = pandas.DataFrame({"#contig":xdf['contig'],"start":xdf['start'],"stop":xdf['stop'],"name":xdf['name']});
	
	return(df_matrix);

def read_result(xinput):
	global gene_list;
	index, path = xinput;
	
	df_ind = pandas.read_csv(path, sep="\t", index_col=False);
	if "bam" not in df_ind.columns or "gw_phased" not in df_ind.columns: return([]);
	
	# change 'sample' to 'sample'
	df_ind.columns = ["sample_id" if x == "bam" else x for x in df_ind.columns];

	out= [];

	for xsample in set(df_ind['sample_id'].tolist()):
		result_all = [];
		result_phased = [];
		df_sample = df_ind[(df_ind.sample_id == xsample)]
		# make sure number and order of genes is the same
		if df_sample['name'].tolist() == gene_list:
			for gw_phased, aCount, bCount in zip(df_sample['gw_phased'].tolist(),df_sample['aCount'].tolist(),df_sample['bCount'].tolist()):
				result_all.append(str(aCount)+"|"+str(bCount));

				if int(gw_phased) == 1:
					result_phased.append(str(aCount)+"|"+str(bCount));
				else:
					result_phased.append("0|0");

			sample = df_sample['sample_id'].tolist()[0];

			result_all = [sample] + result_all;
			list_to_file(result_all, "tmp/all/"+sample+".txt");

			result_phased = [sample] + result_phased;
			list_to_file(result_phased, "tmp/gw_phased/"+sample+".txt");

			out.append(sample);
		else:
			print("ERROR - "+path+":"+xsample+" genes are not in correct order...")

	return(out);

def parallelize(function, threads, pool_input):
	if len(pool_input) > 0:
		threads = min([len(pool_input),threads]);
		if threads > 1:
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

def list_to_file(list,path):
	in_stream = open(path, "w");
	
	for item in list:
		in_stream.write(str(item)+"\n");
	
	in_stream.close();

if __name__ == "__main__":
	main();