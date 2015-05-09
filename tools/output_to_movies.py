import os
from sys import argv

#Ouput_list_name is a text file with the name of an output per line
#order is irrelevant

script, output_list_name, delay = argv

fil = open(str(output_list_name), 'r')

print "-"*40
print "Creating plots"
print "-"*40

for line in fil:
	#Executing from python the output_to_plots.py for each line
	exec_string = "python output_to_plots.py "+str(line)[0:-1]
	print exec_string
	store_for_later = str(line)[0:10]
	os.system(exec_string)
	del exec_string

print "-"*40
print "Creating movies"
print "-"*40


for q in ["d", "p", "u", "v"]:
	exec_string = "convert " + store_for_later+"*.vts_"+q+".png -delay "+str(int(delay)) + " movie_"+q+".mpg"
	print exec_string
	os.system(exec_string)
	del exec_string

fil.close()

print "-"*40
print "Done"
print "-"*40


