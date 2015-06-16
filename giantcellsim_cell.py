import random as rand

class Cell:

	def __init__(self,strands,motif,max_strand_nr,nr_motifs,has_motif,nr_bases):
		self.strands = strands
		self.motif = motif
		self.max_strand_nr = max_strand_nr
		self.nr_motifs = self.motif_count()
		self.has_motif = self.check_for_motif()
		self.nr_bases = self.find_nr_bases()


	def check_for_motif(self):
		if self.nr_motifs > 0:
			return True
		else:
			return False

	def motif_count(self):
		motif_count = 0
		for strand_iterator in range(self.nr_strands()):
			if str(self.motif) in self.strands[strand_iterator]:
				motif_count += 1
		return motif_count

	def update_motifs(self):
		self.nr_motifs = self.motif_count()
		self.has_motif = self.check_for_motif()

	def update_nr_bases(self):
		self.nr_bases = self.find_nr_bases()

	def find_nr_bases(self):
		nr_bases = 0
		for strand in self.strands:
			nr_bases += len(strand)
		return nr_bases

	def nr_strands(self):
		strand_counter = 0
		for strand in self.strands:
			if len(strand) > 0:
				strand_counter += 1
		return strand_counter

	def concentration_bias(self,bias):
		if len(self.strands) > 0:
			motif_concentration = float(self.nr_motifs)/self.nr_strands()
		else:
			motif_concentration = 0
		custom_bias = 0.5+(bias-0.5)*motif_concentration
		return custom_bias


	def grow(self,elong,bias,maxStrandLength):
		strand_counter = 0 
		custom_bias = self.concentration_bias(bias) # call before growth stage in order for motif effect to be consistent
		for strand_iterator in range(self.nr_strands()): 
			uniform_rv = rand.uniform(0,1)
			if uniform_rv < elong and len(self.strands[strand_counter]) < maxStrandLength:
				if self.has_motif == True:
					if rand.uniform(0,1) < custom_bias:
						self.strands[strand_counter] = self.strands[strand_counter] + "0"
					else:
						self.strands[strand_counter] = self.strands[strand_counter] + "1"
				else:
					if rand.uniform(0,1) < 0.5:
						self.strands[strand_counter] = self.strands[strand_counter] + "0"
					else:
						self.strands[strand_counter] = self.strands[strand_counter] + "1"
				strand_counter = strand_counter + 1
			elif uniform_rv < 2*elong: # delete a strand with same probablitity as elongating (or creating new)
				self.strands.pop(strand_counter)


		for empty_iterator in range(self.max_strand_nr-self.nr_strands()):
			if rand.uniform(0,1) < elong:
				if self.has_motif == True:
					if rand.uniform(0,1) < custom_bias:
						self.strands.append("0")
					else:
						self.strands.append("1")
				else:
					if rand.uniform(0,1) < 0.5:
						self.strands.append("0")
					else:
						self.strands.append("1")

		self.update_motifs()
		self.update_nr_bases()


