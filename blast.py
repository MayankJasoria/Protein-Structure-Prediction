from Bio.Blast import NCBIWWW
from Bio.Blast import NCBIXML
from Bio import SeqIO

class BLAST:

    def __init__(self):
        self.E_VALUE_THRESH = 0.04
        self.hitlist_size = 5

    def execute(self, file):
        print("\nRunning BLAST")
        record = SeqIO.read(file, format="fasta")
        result_handle = NCBIWWW.qblast("blastp", "nr", record.format("fasta"), hitlist_size = self.hitlist_size)

        with open("my_blast.xml", "w") as out_handle:
            out_handle.write(result_handle.read())
        result_handle.close()

        result_handle = open("my_blast.xml")

        blast_record = NCBIXML.read(result_handle)

        result = []

        print("COmparing BLAST results with predicted secondary structures")

        for alignment in blast_record.alignments:
            for hsp in alignment.hsps:
                if hsp.expect < self.E_VALUE_THRESH:
                    # print("****Alignment****")
                    one_res = {}
                    one_res["sequence"] = alignment.title
                    one_res["length"] = alignment.length
                    one_res["e value"] = hsp.expect
                    one_res["hsp_query"] = hsp.query[:]
                    one_res["hsp_match"] = hsp.match[:]
                    one_res["hsp_sbjct"] = hsp.sbjct[:]
                    result.append(one_res)

        return str(result)
                    
