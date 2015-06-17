import csv
from giantcellsim_trial import giantcellsim_trial
import itertools
import numpy


def giantcellsim_motifoutput(parameterlist,masterprefix,testprefix,trials,growthIterations,max_strand_nr,maxStrandLength,numCells,numRounds,motif,elong,bias):
	pop_tracker = []

	with open(masterprefix+ testprefix +'_MotifData_motif{motif}_len{maxStrandLength}_bias{bias}_elong{elong}_{trials}trials_numRound{numRounds}.csv'.format(motif = motif, maxStrandLength = maxStrandLength, bias=bias, elong=elong, trials=trials, numRounds=numRounds), 'wb') as f: 
		writer = csv.writer(f)
		writer.writerow(parameterlist)

		for trial in range(trials):
			pop_tracker.append([])
			nr_motifs, nr_strands, nr_cells_with_motif, pop_tracker[trial] = giantcellsim_trial(motif,growthIterations,max_strand_nr,maxStrandLength,numCells,numRounds,elong,bias)

			motif_freq = [motifs / float(total) for motifs,total in itertools.izip(nr_motifs,nr_strands)]
			strands_freq = [strands / float(max_strand_nr*numCells) for strands in nr_strands]
			cells_with_freq = [cells / float(numCells) for cells in nr_cells_with_motif]

			writer.writerow(motif_freq)
			writer.writerow(strands_freq)
			writer.writerow(cells_with_freq)

			if trial == 0:
				motif_freq_aggregate = motif_freq
				strands_freq_aggregate = strands_freq
				cells_with_freq_aggregate = cells_with_freq
				nr_strands_per_time = nr_strands
			else:
				motif_freq_aggregate = [list(round_data) for round_data in zip(motif_freq_aggregate,motif_freq)]
				strands_freq_aggregate = [list(round_data) for round_data in zip(strands_freq_aggregate,strands_freq)]
				cells_with_freq_aggregate = [list(round_data) for round_data in zip(cells_with_freq_aggregate,cells_with_freq)]
				nr_strands_per_time = [list(round_data) for round_data in zip(nr_strands_per_time,nr_strands)]
		
		means = []
		stdevs = [] 

		for iterator in range(3):
			means.append([])
			stdevs.append([])

		for time_point in range(numRounds):
			means[0].append(numpy.mean(numpy.asarray(motif_freq_aggregate[time_point])))
			stdevs[0].append(numpy.std(numpy.asarray(motif_freq_aggregate[time_point])))
			means[1].append(numpy.mean(numpy.asarray(strands_freq_aggregate[time_point])))
			stdevs[1].append(numpy.std(numpy.asarray(strands_freq_aggregate[time_point])))
			means[2].append(numpy.mean(numpy.asarray(cells_with_freq_aggregate[time_point])))
			stdevs[2].append(numpy.std(numpy.asarray(cells_with_freq_aggregate[time_point])))

		for mean_data in means:
			writer.writerow(mean_data)

		for stdev_data in stdevs:
			writer.writerow(stdev_data)
	f.close()

	return pop_tracker, nr_strands_per_time
