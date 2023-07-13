f1 = open('../Trinity/trinity_out_dir.Trinity.fasta', 'r')
f2 = open('quant_out_name.txt', 'r')
f3 = open('quant_out.fasta', 'w')

fasta = {}

for line in f1.readlines():
    if line.startswith(">"):
        name = line[:].rstrip()
        fasta[name] = ''
        continue
    fasta[name] += line.rstrip()

for line in f2.readlines():
    if line[0] == 'T':
        splitline = line.split()
        N = ">" + splitline[0]
        for key in fasta.keys():
            splitline2 = key.split()
            name = splitline2[0]
            if N == name:
                f3.write(key)
                f3.write('\n' + fasta[key] + '\n')

f1.close()
f2.close()
f3.close()
