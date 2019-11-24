from choufasman import ChouFasman
from gor import GOR3
from Z import ZFactor

cf = ChouFasman()
prediction = cf.execute('MKIDAIVGRNSAKDIRTEERARVQLGNVVTAAALHGGIRISDQTTNSVETVVGKGESRVLIGNEYGGKGFWDNHHHHHH')
print('MKIDAIVGRNSAKDIRTEERARVQLGNVVTAAALHGGIRISDQTTNSVETVVGKGESRVLIGNEYGGKGFWDNHHHHHH')
print(prediction)

zf = ZFactor()
zpred = zf.execute('MKIDAIVGRNSAKDIRTEERARVQLGNVVTAAALHGGIRISDQTTNSVETVVGKGESRVLIGNEYGGKGFWDNHHHHHH')
#print('MKIDAIVGRNSAKDIRTEERARVQLGNVVTAAALHGGIRISDQTTNSVETVVGKGESRVLIGNEYGGKGFWDNHHHHHH')
print(zpred)

# BBBBBBBTB AAAAAAAAAAAAAAAAAAAAAAATTBBBBTBTBTBBBBBBBBTB AAAAAATTTTTTAAAATA
tg = GOR3()
tg.load_model()
print(tg.predict('MKIDAIVGRNSAKDIRTEERARVQLGNVVTAAALHGGIRISDQTTNSVETVVGKGESRVLIGNEYGGKGFWDNHHHHHH'))