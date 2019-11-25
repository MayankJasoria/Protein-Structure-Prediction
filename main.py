from choufasman import ChouFasman
from gor import GOR3
from rost import Rost
from Z import ZFactor
from blast import BLAST
from deepGoPlus import DeepGoPlus

class MainApp:

    def __init__(self):
        self.cf = ChouFasman()
        self.zf = ZFactor()
        self.gr = GOR3()
        self.gr.load_model()
        self.rt = Rost()
        self.bt = BLAST()
        self.dgp = DeepGoPlus()

    def predict(self, sequence):

        # similarly add appropriate predictions for rost, blast, dgp

        #>sp|P04637|P53_HUMAN Cellular tumor antigen p53 OS=Homo sapiens OX=9606 GN=TP53 PE=1 SV=4

        rost_prediction = self.rt.execute(sequence)
        bt_prediction = ""

        try:
            fs = open("temp.fasta", "w")
            fs.write(">sp|P04637|P53_HUMAN Cellular tumor antigen p53 OS=Homo sapiens OX=9606 GN=TP53 PE=1 SV=4\n")
            fs.write(sequence)
        except:
            bt_prediction = "Could not perform BLAST due to internal error"
        finally:
            fs.close()

        if bt_prediction == "":
            bt_prediction = self.bt.execute("temp.fasta")
        dgp_prediction = self.dgp.execute(sequence)

        cf_prediction = self.cf.execute(sequence)
        
        gor_prediction = self.gr.predict(sequence)
    
        zf_prediction = self.zf.execute(sequence)

        return cf_prediction, gor_prediction, rost_prediction, zf_prediction, bt_prediction, dgp_prediction