from choufasman import ChouFasman
from gor import GOR3
# from rost import Rost
from Z import ZFactor
# from blast import Blast
# from deepGoPlus import deepGoPlus

class MainApp:

    def __init__(self):
        self.cf = ChouFasman()
        self.zf = ZFactor()
        self.gr = GOR3()
        self.gr.load_model()
        #self.rt = Rost()
        #self.bt = Blast()
        #self.dgp = DeepGoPlus()

    def predict(self, sequence):
        cf_prediction = self.cf.execute(sequence)
        
        gor_prediction = self.gr.predict(sequence)
    
        zf_prediction = self.zf.execute(sequence)

        # similarly add appropriate predictions for rost, blast, dgp

        rost_prediction = ""
        bt_prediction = ""
        dgp_prediction = ""

        return cf_prediction, gor_prediction, rost_prediction, zf_prediction, bt_prediction, dgp_prediction