class ZFactor:

    def __init__(self):
        self.hmap = {'A':0.87,'B':0,'C':1.23,'D':-2.46,'E':-2.53,
                'F':2.68,'G':1.01,'H':0.92,'I':2.16,'J':0,
                'K':2.49,'L':2.29,'M':1.71,'N':0.3,'O':0,
                'P':0.9,'Q':0.3,'R':2.99,'S':0.85,'T':0.95,
                'U':0,'V':1.61,'W':2.96,'X':0,'Y':1.67,'Z':0}
        
    def execute(self, seq):
        seq_len = len(seq)
        tot1 = 0
        tot2 = 0
        tot = 0
        for i in range(seq_len):
            tot = tot + self.hmap[seq[i]]
        tot = (tot/seq_len)
        
        for i in range(seq_len):
            if (seq[i] == 'D' or seq[i] == 'E' or seq[i] == 'G' or seq[i] == 'K' or seq[i] == 'N' or seq[i] == 'R' or seq[i] == 'H'):
                tot1 = tot1 + 1
        
        for i in range(seq_len):
            if (seq[i] == 'F' or seq[i] == 'I' or seq[i] == 'L' or seq[i] == 'M' or seq[i] == 'V' or seq[i] == 'Y'):
                tot2 = tot2 + 1
         
        load = tot1/tot2
        ans = ''
        print(tot1)
        print(tot2)
        
        print("R3 = " , load)
        tr = -0.345*load + 0.6*tot 
        print("Discriminant Factor = " , tr)
        print("H_phi = " , tot)
        
        #print("Is it a membrance protein? [Y/N]")
        #ch = input()
        
        ans = ans + "**If it's a membrane protein: \n"
        if(tr >= 0.38 and tr <= 0.66):
            ans = ans + "It's probably an internal membrane protein."
            
        elif(tr >= -.04 and tr <= 0.3):
            ans = ans + "It's probably an external membrane protein."
        else: ans = ans + "Can't be determined!!!" 
        
        ans = ans + "**If it's a non-membrane protein: \n"
        if(tr>= -.03 and tr<= 0.3):
            ans = ans + "It's probably a non membrane membrane protein."         
        else: ans = ans + "Can't be determined!!!"
        return ans