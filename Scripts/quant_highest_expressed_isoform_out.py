f = open('../Salmon/salmon_quant/quant.sf', 'r')
output = open('quant_out.sf', 'w')

KEY1 = []
KEY2 = []
VALUE = []
K_V = []
DICT = {}
n = 0

for line in f:
    if line[0] == 'T':
        splitline = line.split()
        ID = splitline[0]
        TPM = splitline[3]
        ID_TPM = {ID:TPM}
        DICT.update(ID_TPM)

for key in DICT.keys():
    t1 = key.split('i')
    t2 = t1[0]
    k = key
    KEY1.append(k)
    KEY2.append(t2)
KEY2.append('!')

for value in DICT.values():
    v = float(value)
    VALUE.append(v)

while True:
    if KEY2[n+1] == '!':break
    elif KEY2[n] == KEY2[n+1]:
        if VALUE[n] > VALUE[n+1]:
            del KEY1[n+1]
            del KEY2[n+1]
            del VALUE[n+1]
        else:
            del KEY1[n]
            del KEY2[n]
            del VALUE[n]
    else:
        n = n+1

K_V = [list(i) for i in zip(KEY1,VALUE)]

output.write('Name\tTPM\n')
for row in K_V:
    rowtxt = '{}    {}'.format(row[0],row[1])
    output.write(rowtxt)
    output.write("\n")
output.close()
f.close()

