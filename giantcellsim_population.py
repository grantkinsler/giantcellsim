from giantcellsim_cell import Cell

class Population:

	def __init__(self,cells,nr_motifs,nr_strands,nr_cells_with_motif):
		self.cells = cells
		self.nr_motifs = self.count_motifs()
		self.nr_strands = self.count_strands()
		self.nr_cells_with_motif = self.count_cells_with_motif()

	def populate(self,numCells,motif,max_strand_nr):
		for cell_iterator in range(numCells):
			self.cells.append(Cell([],motif,max_strand_nr,'empty','empty','empty'))

	def update_counters(self):
		self.nr_motifs = self.count_motifs()
		self.nr_strands = self.count_strands()
		self.nr_cells_with_motif = self.count_cells_with_motif()

	def count_motifs(self):
		motif_count = 0
		for cell in self.cells:
			motif_count = motif_count + cell.motif_count()

		return motif_count

	def count_strands(self):
		strand_count = 0
		for cell in self.cells:
			strand_count = strand_count + cell.nr_strands()
		return strand_count

	def count_cells_with_motif(self):
		cell_count = 0
		for cell in self.cells:
			if cell.has_motif == True:
				cell_count += 1

		return cell_count

	def returncontents(self):
		contents = []

		for cell in self.cells:
			contents.append(cell.strands)

		return contents
		