from Bio.Blast import NCBIWWW
from Bio.Blast import NCBIXML
from Bio import SeqIO
record = SeqIO.read(r'C:\Users\Jayesh\Desktop\ins.fasta', format="fasta")
result_handle = NCBIWWW.qblast("blastp", "nr", record.format("fasta"), hitlist_size = 3)

with open("my_blast.xml", "w") as out_handle:
    out_handle.write(result_handle.read())
result_handle.close()

result_handle = open("my_blast.xml")

blast_record = NCBIXML.read(result_handle)

E_VALUE_THRESH = 0.04

for alignment in blast_record.alignments:
     for hsp in alignment.hsps:
         if hsp.expect < E_VALUE_THRESH:
             print("****Alignment****")
             print("sequence:", alignment.title)
             print("length:", alignment.length)
             print("e value:", hsp.expect)
             print(hsp.query[0:75] + "...")
             print(hsp.match[0:75] + "...")
             print(hsp.sbjct[0:75] + "...")
             print(" ")
