#!/usr/bin/env python

# do logos on Pro, Gly, Asp, rest...

f = open('CSSs.txt', 'r')
pro = open('pro.txt', 'w')
gly = open('gly.txt', 'w')
asp = open('asp.txt', 'w')
glu = open('glu.txt', 'w')
rest = open('rest.txt', 'w')
rest_fa = open('rest.fa', 'w')

# Pro CCN
# Gly GGN
# Asp GA[TC]

p2 = 0
gy2 = 0
a2 = 0
gu2 = 0
rest2 = 0


i = 1
for line in f:
	seq = line.strip()
	p = seq[15:18]
	a = seq[18:21]
	min2 = seq[9:11]
	if p[0:2] == "CC":
		pro.write(seq+"\n")
		if min2 == "GA":
			p2 = p2 + 1
	elif p[0:2] == "GG":
		gly.write(seq+"\n")
		if min2 == "GA":
			gy2 = gy2 + 1
	elif p == "GAT" or p == "GAC":
		asp.write(seq+"\n")
		if min2 == "GA":
			a2 = a2 + 1
	elif a == "GAA" or a == "GAG":
		glu.write(seq+"\n")
		if min2 == "GA":
			gu2 = gu2 + 1
	else:
		rest_fa.write(">"+str(i)+"\n")
		rest_fa.write(seq+"\n")
		rest.write(seq+"\n")
		if min2 == "GA":
			rest2 = rest2 + 1
	i = i+1


f.close()
pro.close()
gly.close()
asp.close()
glu.close()
rest.close()
rest_fa.close()

print "Pro and neg.charged at -2:"
print p2
print "Gly and neg.charged at -2:"
print gy2
print "Asp and neg.charged at -2:"
print a2
print "Glu and neg.charged at -2:"
print gu2
print "Rest and neg.charged at -2:"
print rest2






