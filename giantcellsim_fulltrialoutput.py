
import csv

def giantcellsim_fulltrialoutput(parameterlist,masterprefix,testprefix,trial_data,trials,growthIterations,max_strand_nr,maxStrandLength,numCells,numRounds,motif,elong,bias):

	with open(masterprefix + testprefix +'_FullTrial1Data_motif{motif}_len{maxStrandLength}_bias{bias}_elong{elong}_{trials}trials_numRound{numRounds}.csv'.format(motif = motif, maxStrandLength = maxStrandLength, bias= bias, elong=elong, trials=trials, numRounds=numRounds), 'wb') as f:

		strand_writer = csv.writer(f, quotechar="'", quoting=csv.QUOTE_ALL) # making sure we put quotes around our strings so they're not read as numbers
		parameter_writer = csv.writer(f) # we don't need the quotes for the parameter list

		parameter_writer.writerow(parameterlist)

		for time_point in range(numRounds):
			for cell in range(numCells):
				strand_writer.writerow(trial_data[time_point][cell])

	f.close()
