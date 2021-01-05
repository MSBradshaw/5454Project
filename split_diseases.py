num = 0
for line in open('DisGeNET_genesets.txt','r'):
	with open('SubSets/disease_subset_' + str(num) + '.txt','w') as out:
		out.write(line)
	num += 1
