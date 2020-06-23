### get random reads - plot logo

import sys
import gzip
from random import *


fout = open(sys.argv[-1], 'w')
l = sys.argv[-2]

# select reads of length 28
reads = []
for line in sys.stdin:
	if line.split()[5] == l:
		reads.append(line.split()[9])

# pick random fragments
# indices = sample(range(0, len(reads28)), 10000)

# for i in indices:
# 	fout.write(reads28[i]+'\n')

for r in reads:
	fout.write(reads28[r]+'\n')

fout.close()

