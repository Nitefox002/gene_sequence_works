pip install Bio

from Bio import Entrez, SeqIO, Seq

Entrez.email= "124010055@sastra.ac.in"
fh= Entrez.esearch(db= 'nuccore', term= 'NC_001477', retmax= 10)
data= Entrez.read(fh)
fh.close()
#print(data)

idlist= data["IdList"]
handle= Entrez.efetch(db="nuccore", id='NC_001477', rettype="fasta", retmode="text")
records= handle.read()
records= records.split('\n')[1:]
seq= ''.join(records)
print(seq)
#print(type(seq))
#print(len(seq))

from Bio.Seq import Seq
orf= 0
revseq= seq[::-1]
fwd_seq= []
rev_seq= []
for i in range(0,len(seq)):
  dt= seq[i:i+3]
  fwd_seq.append(dt)
  #writefwd= seq[i:i+3]+'\n'
  #temp= i+(i+1)+(i+2)

for i in range(0,len(seq)):
  ds= revseq[i:i+3]
  rev_seq.append(ds)

#Finding no. of Forward and reverse ORFs
print(f"No. of ORFs in Forward and Reverse sequence: {2*len(fwd_seq)}")
print("Forward ORFs: " ,fwd_seq)
print("Reverse ORFs: " ,rev_seq)

#Writing ORFs to files
with open("fwdorf.txt", 'w') as cd:
  for i in fwd_seq:
    cd.write(i+"\n")
    #cd.write("\n")

with open("revorf.txt", 'w') as dd:
  for i in rev_seq:
    dd.write(i+'\n')


#seqobj= Seq(seq)
#print(seqobj)
tc_seq= ""
ts_seq= ""
for i in range(len(fwd_seq)):
  temp1= seq[i:i+3]
  cont= Seq(temp1)
  tranc= cont.transcribe()
  tc_seq= tc_seq + tranc

print("Transcripted Sequence: ",tc_seq)

for i in range(len(tc_seq)):
  temp2= tc_seq[i:i+3]
  cont2= Seq(temp2)
  transl= cont2.translate()
  ts_seq += transl

print("Translated sequence: ", ts_seq)

print("No. of amino acids: ", len(ts_seq))