import sys

if len(sys.argv) != 1+1:
	sys.exit("""
	Usage: %s distribution
	""" % sys.argv[0])

distribution = [list(map(float,i.strip().split())) for i in open(sys.argv[1])]
#          x_mean,       y,  x_min, x_max
x_y = [((i[0]+i[1])/2., i[2], i[0], i[1]) for i in distribution]

for i,j in zip(x_y, sorted(x_y, reverse=True)):
	if i[0]<j[0]:
		print(i[2], i[3], j[1])
	else:
		print(i[2], i[3], i[1])
