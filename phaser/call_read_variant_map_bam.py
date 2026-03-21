import argparse
import sys

import pysam

import read_variant_map


def do_read_variant_map_bam(variant_table, baseq, o, splice, isize_cutoff):
	"""
	Map reads from a BAM stream (stdin) to variant alleles.

	This is a drop-in replacement for the original SAM-text mapper pipeline:
	  samtools view ... | python call_read_variant_map.py

	Here we read BAM directly via pysam to avoid SAM formatting/parsing overhead.
	"""
	args = {
		"variant_table": variant_table,
		"baseq": baseq,
		"o": o,
		"splice": splice,
		"isize_cutoff": isize_cutoff,
	}
	# read_variant_map helper functions (split_read, etc.) rely on a module-level "args" dict.
	read_variant_map.args = args

	stream_out = open(args["o"], "w")
	stream_variants = open(args["variant_table"], "r")

	# Read filtered alignments from stdin (BAM).
	bam_in = pysam.AlignmentFile("-", "rb")
	contigs = list(bam_in.references)
	contig_index = {c: i for i, c in enumerate(contigs)}

	line_variant = read_variant_map.get_next_variant(stream_variants)
	variant_buffer = []

	read_counter = 0
	for read in bam_in.fetch(until_eof=True):
		# Skip headers: pysam yields alignments only.
		if read.is_unmapped:
			continue

		read_chr = read.reference_name
		# pysam is 0-based; phASER mapper expects 1-based POS.
		read_pos = int(read.reference_start) + 1
		template_length = abs(int(read.template_length))

		# Clear variant buffer.
		buffer_remove = []
		for variant_i in range(0, len(variant_buffer)):
			v = variant_buffer[variant_i]
			if v.chr != read_chr:
				# remove variants from previous chromosomes
				if contig_index.get(v.chr, -1) < contig_index.get(read_chr, -1):
					buffer_remove.append(variant_i)
			elif v.pos < read_pos:
				# remove variants behind current read start
				buffer_remove.append(variant_i)
		for index in reversed(buffer_remove):
			del variant_buffer[index]

		if args["isize_cutoff"] != 0 and template_length > args["isize_cutoff"]:
			continue

		# Alignment score (AS tag) is used later for percentile cutoff in phaser.py.
		try:
			alignment_score = int(read.get_tag("AS"))
		except Exception:
			alignment_score = 0

		# Use SAM-style ASCII qualities string (Phred+33) expected by split_read().
		# This avoids per-base Python conversion from query_qualities -> string.
		baseqs = read.qual
		if baseqs is None:
			baseqs = "I" * (len(read.query_sequence) if read.query_sequence else 0)

		bases = read.query_sequence
		cigar = read.cigarstring
		if bases is None or cigar is None:
			continue

		halt = 0
		if line_variant is not None and line_variant.chr != read_chr:
			if line_variant.chr not in contig_index:
				sys.stderr.write(
					"Error, VCF and BAM contigs do not match VCF = %s BAM = %s\n"
					% (line_variant.chr, read_chr)
				)
				sys.exit(1)
			vindex = contig_index[line_variant.chr]
			bindex = contig_index.get(read_chr, -1)
			if vindex < bindex:
				while line_variant is not None and line_variant.chr != read_chr:
					line_variant = read_variant_map.get_next_variant(stream_variants)
			elif bindex > vindex:
				halt = 1

		if halt == 0:
			if line_variant is not None and line_variant.pos < read_pos:
				while (
					line_variant is not None
					and line_variant.pos < read_pos
					and line_variant.chr == read_chr
				):
					line_variant = read_variant_map.get_next_variant(stream_variants)

			alignments = read_variant_map.split_read(read_pos, bases, baseqs, cigar, read.query_name)

			for alignment in alignments:
				# Add variants through the end of this alignment chunk.
				if line_variant is not None:
					alignment_end = (alignment.read_start + read_pos) + len(alignment.pseudo_read)
					while line_variant is not None and line_variant.chr == read_chr and line_variant.pos <= alignment_end:
						variant_buffer.append(line_variant)
						line_variant = read_variant_map.get_next_variant(stream_variants)

				# Only test variants that fall within this alignment chunk. This matches the
				# previous behavior (variants outside return ""), but avoids extra calls.
				alignment_start = alignment.read_start + read_pos
				alignment_end = alignment_start + len(alignment.pseudo_read)
				for xvar in variant_buffer:
					if xvar.pos < alignment_start:
						continue
					if xvar.pos > alignment_end:
						break
					allele = read_variant_map.identify_allele(alignment, read_pos, xvar)
					if allele != "":
						stream_out.write(
							"\t".join(
								[
									read.query_name,
									xvar.id,
									xvar.rs_id,
									allele,
									str(alignment_score),
									xvar.genotype,
									xvar.maf,
								]
							)
							+ "\n"
						)

		read_counter += 1
		if read_counter % 100000 == 0:
			# Match original mapper behavior (stdout). phaser.py redirects stdout to devnull.
			print(
				(
					"               processed %d reads, buffer_size = %d, position = %s:%d"
					% (read_counter, len(variant_buffer), read_chr, read_pos)
				)
			)

	stream_out.close()


def main():
	parser = argparse.ArgumentParser()
	parser.add_argument("--variant_table", type=str, required=True)
	parser.add_argument("--baseq", type=int, default=10)
	parser.add_argument("--o", type=str, required=True)
	parser.add_argument("--splice", type=int, default=1)
	parser.add_argument("--isize_cutoff", type=float, default=0)
	args = parser.parse_args()

	do_read_variant_map_bam(
		variant_table=args.variant_table,
		baseq=args.baseq,
		o=args.o,
		splice=args.splice,
		isize_cutoff=args.isize_cutoff,
	)


if __name__ == "__main__":
	main()
