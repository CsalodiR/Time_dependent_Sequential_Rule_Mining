import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

from modules.outer_func_v5 import *


NoEvents = 10
SeqEndProb = 0.4
StampEndProb = 1.0
NoInPuts = 500
command =["java -jar spmf.jar run TopSeqRules input.txt output.txt 100 10%"][0]
#command =["java -jar spmf.jar run CMRules input_v2.txt output.txt 5% 5%"][0]
#command =["java -jar spmf.jar run CMRules input_v2.txt output.txt 5% 5%"][0]
#command =["java -jar spmf.jar run CMRules input_v3.txt output_v3.txt 5% 5%"][0]
# java -jar spmf.jar run CMRules contextPrefixSpan.txt output.txt 75% 50%


StateTransition = np.random.rand(NoEvents,NoEvents)
StateTransition = StateTransition / np.sum(StateTransition,axis = 0)


InPutSequences = list()



for i in range(NoInPuts):
    
    actseq = np.array([])
    firstEvent = np.floor(np.random.rand(1)*NoEvents)
    # actseq = np.concatenate([actseq, np.array([(i+1)*1000])])
    # actseq = np.concatenate([actseq, np.array([-1])])
    actseq = np.concatenate([actseq, np.array(firstEvent)])
    actseq = np.concatenate([actseq, np.array([-1])])
       
    makeMore = True
    
    j=2
    while makeMore:
        
        makeMoreA = True
        
        while makeMoreA:
            
            if actseq[j-1] == -1:
            
                Probs = StateTransition[:,int(actseq[j-2])]
            else:
                Probs = StateTransition[:,int(actseq[j-1])]
            Probs = np.cumsum(Probs)
        
            moreEventA = np.random.rand(1)
            moreEvent = np.where(Probs > moreEventA)[0][0]
            actseq = np.hstack([actseq, moreEvent])
            
            if np.random.rand(1) < StampEndProb:
                makeMoreA = False
        
        if np.random.rand(1) < SeqEndProb:
            # actseq = np.concatenate([actseq, np.array([-1])])
            # actseq = np.concatenate([actseq, np.array([(i+1)*1000])])
            actseq = np.concatenate([actseq, np.array([-2])])
            makeMore = False
        else:
            actseq = np.concatenate([actseq, np.array([-1])])
        j=j+2    
    InPutSequences.append(actseq)
    


Pattern=list()

for i in range(len(InPutSequences)):
    dum=InPutSequences[i]
    
    
    Pattern.append(np.asarray(dum, dtype = 'int'))
    


writeSPMF3(Pattern)

#%% Exectuting Sequence Mining 

import os

#command =["java -jar spmf.jar run RuleGrowth input.txt output.txt 1% 10%"][0]

os.system(command)

# #%% Read frequent sequences

Rules=[]

Support=[]
Confidence=[]
    
namesOut="output.txt"

text_file = open(namesOut, "r")
lines = text_file.readlines()

for i in range(len(lines)):
    dum=lines[i]
    dum1=re.split(' ==> | #SUP: | #CONF: ',dum)
    
    # Read Antecedent Part
    dum100=dum1[0]
    dum100=dum100.strip().split(',')
    #dum100=dum100.split(',')
    dum107=dum100.copy()
    dum100=list(map(int, dum100))
    
    # Read Consequent Part
    dum101=dum1[1]
    dum101=dum101.strip().split(',')
    #dum101=dum101.split(',')
    dum101=list(map(int, dum101))
    
    dum102=np.block([np.array(dum100),np.array([-1]),np.array(dum101)])
    
    # Read Support
    dum105=dum1[2]
    dum105=dum105.strip().split(',')
    dum105=list(map(int, dum105))
    
    # Read Confidence
    dum106=dum1[3]
    dum106=dum106.strip().split(',')
    dum106=list(map(float, dum106))
    
    
    Rules.append(dum102)
    Support.append(dum105)
    Confidence.append(dum106)
    
#%%

if 1:
    TID_e = list()
    
    for i in range(len(Rules)):
        print(i)
        dum = Rules[i].copy()
        ActTID = list()
        for j in range(len(InPutSequences)):
            #dum1 = InPutCSOut[j][0]
            dum1 = Pattern[j]
    
            #if is_supported2(dum, dum1):
            if is_supported7(dum, dum1):
                ActTID.append(j)
        TID_e.append(ActTID)


#%% Check support Correctness

Selector = 1
i = 1


# Selector = 2
# i = 23

# Selector = 3
# i = 145


dum = Rules[i]
print(dum)
Outer= list()
OuterBNO = list()
patSup = TID_e[i]

indexM = np.where(dum == -1)[0]
indexM = int(indexM)

for k in range(len(patSup)):
    dum1 = InPutSequences[patSup[k]]
    #dum11 = BNOint2str(dum1.tolist(),BNO)
    
    dum11 = list()
    for i in range(len(dum1)):
        dum11.append(int(dum1[i]))
    dum12 = dum11
    #dum11 = np.array(dum11)
    dum13 = ''.join(str(dum12))
    
    Gyujto2 = list()
    Gyujto = list()
    for ki in range(len(dum1)):
        dum3 = dum1[ki]
        
        if Selector == 1:
            if dum3 == dum[0]:
                Gyujto.append(str(dum3))
            elif dum3 == dum[1]:
                Gyujto.append(str("|"))
            elif dum3 == dum[2]:
                Gyujto.append(str(dum3))
            else:
                Gyujto.append(str("_"))
                
        if Selector == 2:
            if dum3 == dum[0]:
                Gyujto.append(str(dum3))
            elif dum3 == dum[1]:
                Gyujto.append(str(dum3))
            elif dum3 == dum[2]:
                Gyujto.append(str("|"))
            elif dum3 == dum[3]:
                Gyujto.append(str(dum3))
            else:
                Gyujto.append(str("_"))
                
        if Selector == 3:
            if dum3 == dum[0]:
                Gyujto.append(str(dum3))
            elif dum3 == dum[2]:
                Gyujto.append(str(dum3))
            elif dum3 == dum[1]:
                Gyujto.append(str("|"))
            elif dum3 == dum[3]:
                Gyujto.append(str(dum3))
            else:
                Gyujto.append(str("_"))


            
    Gyujto= ', '.join(Gyujto)
    OuterBNO.append(dum13)
    Outer.append(Gyujto)





