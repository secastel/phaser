import argparse;
import read_variant_map;

# usage
# samtools view -q 255 -L variants.bed NA06986.2.M_111215_4.bam | python2.7 call_read_variant_map.py

def main():
	#Arguments passed 
	parser = argparse.ArgumentParser()
	# required
	parser.add_argument("--variant_table", type=str, required=True)
	parser.add_argument("--baseq", type=int, default=10)
	parser.add_argument("--o", type=str, required=True)
	parser.add_argument("--splice", type=int, default=1)
	parser.add_argument("--isize_cutoff", type=float, default=0)
	
	args = parser.parse_args()
	
	read_variant_map.do_read_variant_map(args.variant_table, args.baseq, args.o, args.splice, args.isize_cutoff);
	
if __name__ == "__main__":
	main();
