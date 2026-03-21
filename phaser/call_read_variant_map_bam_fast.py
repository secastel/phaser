import argparse
import pysam


class VariantRec:
	__slots__ = ("pos", "id", "rs_id", "genotype", "maf")

	def __init__(self, pos, vid, rs_id, genotype, maf):
		self.pos = pos
		self.id = vid
		self.rs_id = rs_id
		self.genotype = genotype
		self.maf = maf


def load_variant_table(variant_table_path):
	"""
	Load a per-chromosome phASER variant table into a list sorted by 1-based POS.

	Expected columns (tab-delimited):
	  chr  pos  unique_id  rs_id  alleles  ref_len  genotype  maf
	"""
	variants = []
	with open(variant_table_path, "r") as f:
		for line in f:
			fields = line.rstrip("\n").split("\t")
			if len(fields) < 8:
				continue
			pos = int(fields[1])
			variants.append(
				VariantRec(
					pos=pos,
					vid=fields[2],
					rs_id=fields[3],
					genotype=fields[6],
					maf=fields[7],
				)
			)
	variants.sort(key=lambda v: v.pos)
	return variants


def do_read_variant_map_bam_fast(variant_table, baseq, o, isize_cutoff):
	"""
	Fast SNP-only mapper for phASER.

	This avoids the older split_read()/identify_allele() path and instead uses
	pysam's CIGAR-aware reference-position mapping to extract alleles.

	Note: phASER's default behavior excludes indels from the variant table, so
	alleles are single bases and this mapper is correct for the default settings.
	"""
	variants = load_variant_table(variant_table)
	var_i = 0
	nvars = len(variants)

	stream_out = open(o, "w")

	# Read filtered alignments from stdin (BAM).
	bam_in = pysam.AlignmentFile("-", "rb")

	read_counter = 0
	for read in bam_in.fetch(until_eof=True):
		if read.is_unmapped:
			continue

		template_length = abs(int(read.template_length))
		if isize_cutoff != 0 and template_length > isize_cutoff:
			continue

		bases = read.query_sequence
		if bases is None:
			continue

		# Base qualities: array('B') of Phred scores or None.
		quals = read.query_qualities

		# Alignment score (AS tag) is used later for percentile cutoff in phaser.py.
		try:
			alignment_score = int(read.get_tag("AS"))
		except Exception:
			alignment_score = 0

		# Advance global variant index: reads are coordinate-sorted, so once a read start
		# has passed a variant, no future read can cover it.
		read_start = int(read.reference_start) + 1  # 1-based
		while var_i < nvars and variants[var_i].pos < read_start:
			var_i += 1

		# Parse the CIGAR once to build:
		# 1) match segments: reference intervals that map 1:1 to query positions
		# 2) insertions anchored to the previous reference base (so SNP allele can be "CT", etc.)
		match_segments = []
		insert_after = {}  # key: 1-based reference pos of anchor base; val: inserted seq (masked)

		rpos1 = read_start  # 1-based next reference position to consume
		qpos0 = 0  # 0-based next query position to consume
		seen_ref_skip = False
		for op, length in (read.cigartuples or []):
			# Operations: https://samtools.github.io/hts-specs/SAMv1.pdf
			# 0 M, 1 I, 2 D, 3 N, 4 S, 5 H, 6 P, 7 =, 8 X
			if op in (0, 7, 8):  # consumes query and reference (aligned bases)
				match_segments.append((rpos1, rpos1 + length - 1, qpos0))
				rpos1 += length
				qpos0 += length
			elif op == 1:  # insertion to the reference (consumes query only)
				# Match legacy mapper behavior: insertions that occur after a ref-skip (N) are
				# not appended to alleles due to how split_read() resets coordinates.
				# This preserves output compatibility for spliced reads.
				if not seen_ref_skip:
					anchor = rpos1 - 1
					if anchor >= 1:
						ins_seq = bases[qpos0 : qpos0 + length]
						if ins_seq:
							if quals is not None and baseq > 0:
								# Mask low-quality inserted bases to 'N', matching legacy behavior.
								ins_chars = []
								for i, ch in enumerate(ins_seq):
									ins_chars.append(ch if quals[qpos0 + i] >= baseq else "N")
								ins_seq = "".join(ins_chars)
							insert_after[anchor] = insert_after.get(anchor, "") + ins_seq
				qpos0 += length
			elif op == 2:  # deletion consumes reference only
				rpos1 += length
			elif op == 3:  # ref-skip (splice) consumes reference only
				seen_ref_skip = True
				rpos1 += length
			elif op == 4:  # soft clip consumes query only
				qpos0 += length
			else:
				# H/P and other ops consume neither or are not relevant here.
				pass

		# Emit alleles by iterating only over variants within aligned match segments.
		vi = var_i
		for ref_start, ref_end, qstart in match_segments:
			while vi < nvars and variants[vi].pos < ref_start:
				vi += 1
			while vi < nvars and variants[vi].pos <= ref_end:
				v = variants[vi]
				qidx = qstart + (v.pos - ref_start)
				if qidx < 0 or qidx >= len(bases):
					vi += 1
					continue

				base = bases[qidx]
				if quals is not None and quals[qidx] < baseq:
					base = "N"

				allele = base + insert_after.get(v.pos, "")
				if allele != "N":
					stream_out.write(
						"\t".join(
							[
								read.query_name,
								v.id,
								v.rs_id,
								allele,
								str(alignment_score),
								v.genotype,
								v.maf,
							]
						)
						+ "\n"
					)

				vi += 1

		read_counter += 1
		if read_counter % 100000 == 0:
			# Match original mapper behavior (stdout). phaser.py redirects stdout to devnull.
			print("               processed %d reads" % read_counter)

	stream_out.close()


def main():
	parser = argparse.ArgumentParser()
	parser.add_argument("--variant_table", type=str, required=True)
	parser.add_argument("--baseq", type=int, default=10)
	parser.add_argument("--o", type=str, required=True)
	parser.add_argument("--isize_cutoff", type=float, default=0)
	args = parser.parse_args()

	do_read_variant_map_bam_fast(
		variant_table=args.variant_table,
		baseq=args.baseq,
		o=args.o,
		isize_cutoff=args.isize_cutoff,
	)


if __name__ == "__main__":
	main()
