from giantcellsim_cell import Cell
from giantcellsim_population import Population
import random as rand
from copy import deepcopy
from copy import copy

def giantcellsim_trial(motif,growthIterations,max_strand_nr,maxStrandLength,numCells,numRounds,elong,bias):

	population = Population([],'empty','empty','empty')

	population.populate(numCells,motif,max_strand_nr)

	# counter lists
	nr_motifs = []
	nr_strands_used = []
	nr_cells_with_motif = []
	population_tracker = []

	for time in range(numRounds):
		for growth_iter in range(growthIterations):
			for cell_iterator in range(numCells):
				population.cells[cell_iterator].grow(elong,bias,maxStrandLength)

		population.update_counters()

		nr_motifs.append(copy(population.nr_motifs))
		nr_strands_used.append(copy(population.nr_strands))
		nr_cells_with_motif.append(copy(population.nr_cells_with_motif))
		population_tracker.append(deepcopy(population.returncontents()))

	return nr_motifs, nr_strands_used, nr_cells_with_motif, population_tracker