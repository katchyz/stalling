## split disorder into tx files


txs = []
f_tx = open("tx.txt", "r")
for line in f_tx:
	tx = line.strip()
	txs.append(tx)


f_tx.close()


f_dis = open("css_scores.disembl", "r")
counter = -1
for line in f_dis:
	if line.startswith("#"):
		# save previous
		if counter >= 0:
			fout = open(txs[counter]+".txt", "w")
			fout.write(scores)
			fout.close()
		counter = counter + 1
		scores = ''
	else:
		l =  str(line.split()[1])+"\t"+str(line.split()[2])+"\t"+str(line.split()[3])+"\n"
		scores = scores + l




fout = open(txs[counter], "w")
fout.write(str(scores))
fout.close()

f_dis.close()