import sys

if len(sys.argv) != 1+2:
	sys.exit("""
	Usage: %s distribution1 distribution2
	Do 'distribution1 - distribution2' :)
	""" % sys.argv[0])

d1 = [list(map(float,i.strip().split())) for i in open(sys.argv[1])]
d2 = [list(map(float,i.strip().split())) for i in open(sys.argv[2])]

#          x_mean,       y,  x_min, x_max
x_y1 = [((i[0]+i[1])/2., i[2], i[0], i[1]) for i in d1]
x_y2 = [((i[0]+i[1])/2., i[2], i[0], i[1]) for i in d2]

max1 = max([i for i in x_y1[1]])
max2 = max([i for i in x_y2[1]])
ratio = max1/max2

if len(x_y1) != len(x_y2):
	sys.exit("ERROR: distributions have different sizes!")

for n,i in enumerate(x_y1):
	x1 = i[0]
	if x_y2[n][0] != x1:
		sys.exit("ERROR: bin not found in distribution2: %s-%s" % (i[2], i[3]))
	diff = i[1]-x_y2[n][1]
	if diff > 0:
		print(i[2], i[3], diff)
	else:
		print(i[2], i[3], 0)
