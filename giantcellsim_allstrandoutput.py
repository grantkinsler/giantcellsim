
import csv
import collections
import numpy
import itertools

def makeKeyorder(maxStrandLength,motif): # create the keyorder for the order of our dictionary
	
	keyorder = []

	motifcounter = 0

	for n in range(maxStrandLength): # create order of keys for the ordered dictionary
		for key in itertools.product(range(2),repeat = n+1):
			mod_key = str(key).strip(" ,(),','").replace(", ", "")
			if motif in mod_key:
				keyorder.insert(motifcounter,mod_key) # place all motif strands first in the order
				motifcounter = motifcounter + 1
			else:
				keyorder.append(mod_key) # place the rest in the order they are created

	return keyorder

def giantcellsim_allstrandoutput(parameterlist,masterprefix,testprefix,pop_tracker,nr_strands_per_time,trials,growthIterations,max_strand_nr,maxStrandLength,numCells,numRounds,motif,elong,bias):

	keyorder = makeKeyorder(maxStrandLength,motif)

	time_trial_dict = []
	dict_per_time = []
	for trial in range(trials):
		time_trial_dict.append([])
		for time_point in range(numRounds):
			temp_dict = {}
			if trial == 0:
				temp2_dict = {}
			for key in keyorder:
				temp_dict[key] = 0
				if trial == 0:
					temp2_dict[key] = []
			if trial == 0:
				dict_per_time.append(collections.OrderedDict(sorted(temp2_dict.items(), key = lambda i:keyorder.index(i[0]))))
			time_trial_dict[trial].append(collections.OrderedDict(sorted(temp_dict.items(), key = lambda i:keyorder.index(i[0]))))

	for trial in range(trials):
		for time_point in range(numRounds):
			for cell in range(len(pop_tracker[trial][time_point])):
				for strand in pop_tracker[trial][time_point][cell]:
					time_trial_dict[trial][time_point][strand] = time_trial_dict[trial][time_point][strand] + 1

	for trial in range(trials):
		for time_point in range(numRounds):
			for key, value in time_trial_dict[trial][time_point].iteritems():
				dict_per_time[time_point][key].append(int(value)/float(nr_strands_per_time[time_point][trial]))

	stdev_dict = [] 
	mean_dict = []

	for time_point in range(numRounds):
		temp_mean = {}
		temp_stdev = {}
		for key in keyorder:
			temp_mean[key] = numpy.mean(dict_per_time[time_point][key])
			temp_stdev[key] = numpy.std(dict_per_time[time_point][key])
		mean_dict.append(collections.OrderedDict(sorted(temp_mean.items(), key = lambda i:keyorder.index(i[0]))))
		stdev_dict.append(collections.OrderedDict(sorted(temp_stdev.items(), key = lambda i:keyorder.index(i[0]))))


	with open(masterprefix + testprefix +'_AllStrandData_motif{motif}_len{maxStrandLength}_bias{bias}_elong{elong}_{trials}trials_numRound{numRounds}.csv'.format(motif = motif, maxStrandLength = maxStrandLength, bias= bias, elong=elong, trials=trials, numRounds=numRounds), 'wb') as f:
		parameter_writer = csv.writer(f)
		dict_writer = csv.DictWriter(f,mean_dict[0].keys())
		dict_header_writer = csv.DictWriter(f,mean_dict[0].keys(),quotechar="'", quoting=csv.QUOTE_ALL)

		parameter_writer.writerow(parameterlist)

		dict_header_writer.writeheader()
		for time_point in range(numRounds):
			dict_writer.writerow(mean_dict[time_point])
		for time_point in range(numRounds):
			dict_writer.writerow(stdev_dict[time_point])

	f.close()
