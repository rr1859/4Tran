module load python3
module loada biopython

python3
from Bio import SeqIO;import re;
import numpy as np;
from Bio.Seq import Seq
import sys, getopt
re1=''
re2=''
filename=''
myopts, args = getopt.getopt(sys.argv[1:], "i:1:2:")

for o, a in myopts:
    if o == '-i':
        filename=a
    elif o == '-1':
        re1=a
    elif o == '-2':
        re2=a

record = SeqIO.read(open(filename), "fasta");
record = SeqIO.read("IAPEz-int_consensus.fasta", "fasta");
re1_re1_250 = re.findall(r'(?=('+re1+'.(?:.){250,}?'+re1+'))', str(record.seq))
pos_matches = [m.start() for m in re.finditer(r'(?=('+re1+'.(?:.){250,}?'+re1+'))', str(record.seq))];
primary_pos = [m.start() for m in re.finditer(re1, str(record.seq))];
secondary_pos = [m.start() for m in re.finditer(re2, str(record.seq))];
if(primary_pos[0]-secondary_pos[0] > 150 and primary_pos[0] > 250):
    print("beg")
    re1_re1_250.append(str(record.seq)[secondary_pos[0]:primary_pos[0]])
    pos_matches.append(1)
elif(secondary_pos[-1]-primary_pos[-1] > 150 and len(str(record.seq))-primary_pos[-1] > 250):
    print("end")
    re1_re1_250.append(str(record.seq)[primary_pos[-1]:secondary_pos[-1]+len(re1)])
    pos_matches.append(len(str(record.seq)))

valid_primer_seq = []
valid_primer_seq_pos = []
valid_primer_seq_size = []
for i in range(0, len(re1_re1_250)):
    s = re1_re1_250[i]
    s_split_re2 = re.split(re2,s)
    s_split_re1 = re.split(re1,s)
    if(len(s_split_re2) > 1 and len(s_split_re1) >= 2):
        print("yes")
        for s_2 in [s_split_re2[0]+re2, re2+s_split_re2[-1]]:
            #print('primer '+s_2)
            if(len(s_2) > 150):
                #print('valid '+str(len(s_2))+ ' '+s_2)
                valid_primer_seq.append(s_2)
                valid_primer_seq_pos.append(pos_matches[i])
                valid_primer_seq_size.append(len(s))


re1_onleft = np.matrix(['','','','','','',''])
re1_onright = np.matrix(['','','','','','',''])

for i in range(0, len(valid_primer_seq)):
    v_s = valid_primer_seq[i]
    if(v_s.find(re1) == 0):
        left_seq = v_s[-100:]+v_s[0:20]
        rev_comp_rightprimer = Seq(left_seq[-20:])
        re1_onleft = np.vstack((re1_onleft, ['left'+str(i),left_seq,str(rev_comp_rightprimer.reverse_complement()),str(valid_primer_seq_size[i]),v_s, str(len(v_s)),str(valid_primer_seq_pos[i])]))
    else:
        right_seq = v_s[-20:]+v_s[0:100]
        left_primer = right_seq[0:20]
        re1_onright = np.vstack((re1_onright, ['right'+str(i),right_seq,left_primer,str(valid_primer_seq_size[i]),v_s, str(len(v_s)),str(valid_primer_seq_pos[i])]))

np.savetxt('re1_onright.txt', re1_onright[1:,:], delimiter='\t', fmt='%s')
np.savetxt('re1_onleft.txt', re1_onleft[1:,:], delimiter='\t', fmt='%s')

print('\n--files re1_onright.txt and re1_onleft.txt saved\n')
