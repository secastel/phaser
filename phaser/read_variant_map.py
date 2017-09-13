import sys;

def do_read_variant_map(variant_table, baseq, o, splice, isize_cutoff):
	global args;
	args = {};
	args['variant_table'] = variant_table;
	args['baseq'] = baseq;
	args['o'] = o;
	args['splice'] = splice;
	args['isize_cutoff'] = isize_cutoff;
	
	stream_out = open(args['o'], "w");
	stream_variants = open(args['variant_table'],"r");
	#[output['chr'],str(output['pos']-1),str(output['pos']),output['id']

	contigs = [];
	
	line_variant = get_next_variant(stream_variants);

	variant_buffer = [];
	as_column = -1;

	read_counter = 0;
	
	for line in sys.stdin:
		# first read the order of chromosomes
		read_columns = line.rstrip().split("\t");
		if read_columns[0][0:3] == "@SQ":
			contigs.append(line.rstrip().split("\t")[1].split(":")[1]);
		elif read_columns[0][0:1] != "@":
			#ERR188213.33037635	355	1	11673	3	76M	=	11784	352	TGAAGCCCTGGAGATTCTTATTAGTGATTTGGGCTGGGGCCTGGCCATGTGTATTTTTTTAAATTTCCACTGATGA	C@CAFFFFHHHDHFHFE>HHHIJICHIGIGGGIEEHIEGGFFHGHIGFGFAB9BFHGIIEFFDFFFFFDCCEEEEA	PG:Z:MarkDuplicates	RG:Z:id	NH:i:2	HI:i:2	nM:i:0	AS:i:152	
			halt = 0;
			read_chr = read_columns[2];
			read_pos = int(read_columns[3]);
			template_length = abs(int(read_columns[8]));
			
			# clear variant buffer
			buffer_remove = [];
			for variant_i in range(0, len(variant_buffer)):
				if variant_buffer[variant_i].chr != read_chr:
					# remove variants from buffer that are from previous chromosomes
					if contigs.index(variant_buffer[variant_i].chr) < contigs.index(read_chr):
						buffer_remove.append(variant_i);
				elif variant_buffer[variant_i].pos < read_pos:
					# remove variants from buffer that are behind the minimum read position on the current chromosome
					buffer_remove.append(variant_i);
			
			for index in reversed(buffer_remove):
				del variant_buffer[index];
			
			if (args['isize_cutoff'] == 0 or template_length <= args['isize_cutoff']):
				# get the AS
				alignment_score = "";
		
				if as_column == -1:
					for i in range(11,len(read_columns)):
						if read_columns[i].startswith("AS:"):
							alignment_score = int(read_columns[i].split(":")[2]);
							as_column = i;
				else:
					alignment_score = int(read_columns[as_column].split(":")[2]);
				
				#Set as_column back to initial value, in case not all reads have equal number fields
				as_column = -1;
		
				# BAM and variant not on same chromosome
				if line_variant != None:
					if line_variant.chr != read_chr:
						if line_variant.chr not in contigs:
							print("Error, VCF and BAM contigs do not match VCF = %s BAM = %s"%(line_variant.chr,read_chr));
							sys.exit(1)
						else:
							vindex = contigs.index(line_variant.chr);
							bindex = contigs.index(read_chr);
				
							if vindex < bindex:
								# need to seek in variants
								while line_variant.chr != read_chr:
									line_variant = get_next_variant(stream_variants);
									if line_variant == None:
										break;
							elif bindex > vindex:
								# need to seek bam
								halt = 1;
		
				# variant file and BAM are on the same chromosome
				if halt == 0:
					if line_variant != None:
						if line_variant.pos < read_pos:
							# seek the variants
							while line_variant.pos < read_pos and line_variant.chr == read_chr:
								line_variant = get_next_variant(stream_variants);
								if line_variant == None:
									break;
					
					alignments = split_read(read_pos, read_columns[9],read_columns[10],read_columns[5], read_columns[0]);
					
					for alignment in alignments:
						
						intersected_variants = []
						for variant_i in range(0, len(variant_buffer)):
							if variant_buffer[variant_i].pos >= (alignment.read_start + read_pos) and variant_buffer[variant_i].pos <= (alignment.read_start + read_pos)+len(alignment.pseudo_read):
								# otherwise record an intersection
								intersected_variants.append(variant_buffer[variant_i]);
					
						if line_variant != None:
							while line_variant.pos <= (alignment.read_start + read_pos)+len(alignment.pseudo_read) and line_variant.chr == read_chr:
								intersected_variants.append(line_variant);
								variant_buffer.append(line_variant);
								line_variant = get_next_variant(stream_variants);
								if line_variant == None:
									break;
				
						for xvar in variant_buffer:
							allele = identify_allele(alignment, read_pos, xvar);
							if allele != "":
								stream_out.write("\t".join([read_columns[0],xvar.id,xvar.rs_id,allele,str(alignment_score),xvar.genotype,xvar.maf])+"\n");
					
					
		read_counter += 1;
		
		if read_counter%100000 == 0:
			print("               processed %d reads, buffer_size = %d, position = %s:%d"%(read_counter,len(variant_buffer),read_chr,read_pos));
	stream_out.close();
	
class variant:
	def __init__(self, variant_columns):
		self.chr = variant_columns[0]
		self.pos = int(variant_columns[1]);
		self.id = variant_columns[2];
		self.rs_id = variant_columns[3];
		self.alleles = variant_columns[4].split(",");
		self.ref_length = int(variant_columns[5]);
		self.genotype = variant_columns[6];
		self.maf = variant_columns[7];
	
	def __str__(self):
		return("\t".join(map(str,[self.chr,self.pos,self.id])));
				
class alignment:
	def __init__(self, input):
		self.read_start = input[0]
		self.read_stop = input[1];
		self.pseudo_read = input[2];
		self.insertions = input[3];
		self.read_pos = input[4];
		self.read_id = input[5];
		
	def genome_start(self):
		return(self.read_pos + self.read_start);
	
	def __str__(self):
		return("\t".join(map(str,[self.read_start,self.read_stop,self.pseudo_read,self.insertions,self.read_pos])));

def get_next_variant(stream_variants):
	buffer = stream_variants.readline();
	if buffer != "":
		xvar = variant(buffer.rstrip().split("\t"));
	else:
		xvar = None;
		#print("               completed mapping reads...")
	
	return(xvar);

def split_read(alignment_pos, bases,baseqs,cigar, read_id):
	global args;
	
	alignments = [];
	
	if args['splice'] == 1 or "N" not in cigar:
		#first filter read by BASEQ
		number_build = "";
		read_seq = "";
		read_pos = 0;
	
		genome_start = 0;
		genome_pos = 0;
	
		for base, quality in zip(bases, baseqs):
			baseq = ord(quality) - 33;
			if baseq >= args['baseq']:
				read_seq += base;
			else:
				read_seq += "N";
		
		bases = read_seq;
		
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
				if c == "M" or c == "X" or c == "=":
					#alignment
					pseudo_read += bases[read_pos:read_pos+seq_len];
					read_pos += seq_len;
					aligned_bases += seq_len;
					genome_pos += seq_len;
				
				elif c == "N":
					#gap
					alignments.append(alignment([genome_start,genome_pos,pseudo_read,insertions,alignment_pos,read_id]));
					genome_pos += seq_len;
					genome_start = genome_pos;
					pseudo_read = "";
					insertions = {};
					
				elif c == "D":
					#deletion
					pseudo_read += "D" * seq_len;
					genome_pos += seq_len;
					
				elif c == "I":
					#insertion
					insertions[genome_pos-1] = bases[read_pos:read_pos+seq_len];
					read_pos += seq_len;
					
				elif c == "S":
					#clipped, not used in alignment
					read_pos += seq_len;

				elif c == "H":
					#hard-clipped, advance neither on read nor reference
					pass

				number_build = "";
		alignments.append(alignment([genome_start,genome_pos,pseudo_read,insertions,alignment_pos,read_id]));
	
	return(alignments);

def identify_allele(alignment, read_coord, xvar):
	map_start = alignment.genome_start();
	
	read_start = (xvar.pos)-map_start;
	read_end = (xvar.pos+xvar.ref_length)-map_start;
	
	if read_start >= 0 and read_end <= len(alignment.pseudo_read):
		read_seq = alignment.pseudo_read[read_start:read_end]
		
		offset = 0;
		for genome_pos in range(read_start,read_end):
			#print(genome_pos);
			if genome_pos in alignment.insertions:
				insert_pos = (genome_pos - read_start) + offset + 1;
				read_seq = read_seq[0:insert_pos] + alignment.insertions[genome_pos] + read_seq[insert_pos:len(read_seq)];
				offset += len(alignment.insertions[genome_pos]);
		
		#remove deletion placeholder
		read_seq = read_seq.replace("D","");
		if read_seq != "N":
			return(read_seq);
	
	return("");
