import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from datetime import date
from datetime import datetime
import copy
import re




def writeSPMF(Pattern):
    names=["input.txt"][0]
    file = open(names,'w')
    
    for i in range(len(Pattern)):
    
        dum=Pattern[i]
    
        for j in range(len(dum)-1):
            file.write(str(dum[j]))
            file.write(' ')
            
        #actseq = np.concatenate([actseq, np.array([-1])])
        #actseq = np.concatenate([actseq, np.array([(i+1)*1000])])
        file.write(str(-1))
        file.write(" ")
        file.write(str((i+1)*1000))
        file.write(" ")
        
        if i == len(Pattern):
            file.write(str(dum[len(dum)-1]))
        else:
            file.write(str(dum[len(dum)-1])+"\n")
        
    file.close()


def readSeqRules():
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
        
        
        Rules.append(dum102) # X => Y <==> [X -1 Y]
        Support.append(dum105)
        Confidence.append(dum106)
        
    return Rules, Support, Confidence


        
def is_supported7(rule, sequence):

    index=int(np.where(rule == -1)[0])
    antecedent = rule[:index]  # Extract the antecedent from the rule
    consequent = rule[index+1::]  # Extract the consequent from the rule
    
    index1=np.in1d(antecedent,sequence)
    index2=np.in1d(consequent,sequence)
    index3=np.all(index1)
    index4=np.all(index2)
    
    if np.logical_and(index3,index4):
        
        indexZ = list()
        for j in range(len(antecedent)):
            dum11 = np.where(sequence == antecedent[j])[0]
            indexZ.append(int(dum11[0]))
            
        indexY = list()
        for j in range(len(consequent)):
            dum11 = np.where(sequence == consequent[j])[0]
            indexY.append(int(dum11[-1]))
       
        if np.max(np.array(indexZ)) < np.min(np.array(indexY)):
            
            isant = True
            indexH=np.ones(len(rule))*(-1)
            for j in range(len(rule)):
                if rule[j] == -1:
                    isant = False
                    continue
                if isant:
                    indexH[j] = (np.where(sequence == rule[j])[0])[0]
                else:
                    indexH[j] = (np.where(sequence == rule[j])[0])[-1]
                
            indexT=np.where(indexH == -1)[0]
            indexH = indexH[np.array([int(indexT-1),int(indexT+1)])]
            
            indexA=np.where(sequence == -1)[0]
            indexA = np.block([indexA,indexH])
            indexA = np.sort(indexA)
            
            #indexZ = list()
            indexZ = np.array([0,0])
            for j in range(len(indexH)):
                indexZ[j]=(np.where(indexA==indexH[j])[0])
            
            indexZ = np.array(indexZ)
            indexZ = np.diff(indexZ)
            
            if indexZ>1:
                
                return True
            else:
                
                return False
            
            
        else:
            return False
            
    else:
        return False
        
        
        
def findTID(Rules, Pattern):
    TID_e = list()
    
    for i in range(len(Rules)):
        
        dum = Rules[i].copy()
        ActTID = list()
        for j in range(len(Pattern)):
            #dum1 = InPutCSOut[j][0]
            dum1 = Pattern[j]
    
            #if is_supported2(dum, dum1):
            if is_supported7(dum, dum1):
                ActTID.append(j)
        TID_e.append(ActTID)
    return TID_e

def separateRule(Rules):
    Antecedent = list()
    Consequent = list()


    for i in range(len(Rules)):
        dum = Rules[i]
        indexA = int(np.where(dum == -1)[0])
        dumAnt = dum[0:indexA]
        dumCon = dum[indexA+1::]
        dumAntLen = len(dumAnt)
        
        Antecedent.append(dumAnt)
        Consequent.append(dumCon)
        
    return Antecedent, Consequent


def calculateDistributionFunctions(Antecedent,Consequent,Pattern,InPutTimeStamps,TID_e):
    from lifelines import KaplanMeierFitter
    kmf = KaplanMeierFitter()

    SurTime = list()
    for i in range(len(Antecedent)):
        
        dumAnt = Antecedent[i].copy()
        dumCon = Consequent[i].copy()
        dumTID = TID_e[i].copy()
        
        dumSurTime = list()
        for j in range(len(dumTID)):
            dumInputS = Pattern[int(dumTID[j])]
            dumInputT = InPutTimeStamps[int(dumTID[j])]
            
            
            dumAntRed = list()
            for k in range(len(dumAnt)):
                dumAntRed.append(np.where(dumInputS == dumAnt[k])[0][0])
            dumAntRed = np.array(dumAntRed)
            dumAntRed = np.max(dumAntRed)
            
            
            dumConRed = list()
            for k in range(len(dumCon)):
                dumConRedA = np.where(dumInputS == dumCon[k])[0]
                moreCheck = True
                m = 0
                while moreCheck:
                    
                    if dumConRedA[m] > dumAntRed:
                        moreCheck = False
                    m = m + 1
                dumConRed.append(dumConRedA[m-1])
                
                
            dumConRed = np.array(dumConRed)
            dumConRed = np.max(dumConRed)
            
            timeT = dumInputT[int(dumConRed)] - dumInputT[int(dumAntRed)]
            timeT = timeT / 365.25
            dumSurTime.append(timeT)
        dumSurTime=np.array(dumSurTime)
        dumSurTime[dumSurTime<0] = 0
        
        
        kmf.fit(dumSurTime,np.ones(len(dumSurTime)))
        x = np.array(kmf.timeline)
        y = np.array(kmf.survival_function_.iloc[:, 0])  # Assuming a single column for the survival function
        
        SurT = list()
        SurT.append(x)
        SurT.append(y)
        SurTime.append(SurT)
        
    return SurTime

def calculateConfidenceFunctions(SurvivalFunctions,Confidence):
    ConfidenceFunctions = list()
    for i in range(len(SurvivalFunctions)):
        dumSurX = SurvivalFunctions[i][0]
        dumSurY = SurvivalFunctions[i][1]
        dumConfStat = Confidence[i]
        dumConfY = (1-dumSurY)*dumConfStat
        
        loader = list()
        loader.append(dumSurX)
        loader.append(dumConfY)
        ConfidenceFunctions.append(loader)
        
    return ConfidenceFunctions

def calculateDistributionFunction(Antecedent,Consequent,Pattern,InPutTimeStamps,TID_e):
    from lifelines import KaplanMeierFitter
    kmf = KaplanMeierFitter()

    SurTime = list()
    for i in range(len(Antecedent)):
        
        dumAnt = Antecedent.copy()
        dumCon = Consequent.copy()
        dumTID = TID_e.copy()
        
        dumSurTime = list()
        for j in range(len(dumTID)):
            dumInputS = Pattern[int(dumTID[j])]
            dumInputT = InPutTimeStamps[int(dumTID[j])]
            
            
            dumAntRed = list()
            for k in range(len(dumAnt)):
                dumAntRed.append(np.where(dumInputS == dumAnt[k])[0][0])
            dumAntRed = np.array(dumAntRed)
            dumAntRed = np.max(dumAntRed)
            
            
            dumConRed = list()
            for k in range(len(dumCon)):
                dumConRedA = np.where(dumInputS == dumCon[k])[0]
                moreCheck = True
                m = 0
                while moreCheck:
                    
                    if dumConRedA[m] > dumAntRed:
                        moreCheck = False
                    m = m + 1
                dumConRed.append(dumConRedA[m-1])
                
                
            dumConRed = np.array(dumConRed)
            dumConRed = np.max(dumConRed)
            
            timeT = dumInputT[int(dumConRed)] - dumInputT[int(dumAntRed)]
            timeT = timeT / 365.25
            dumSurTime.append(timeT)
        dumSurTime=np.array(dumSurTime)
        dumSurTime[dumSurTime<0] = 0
        
        
        kmf.fit(dumSurTime,np.ones(len(dumSurTime)))
        x = np.array(kmf.timeline)
        y = np.array(kmf.survival_function_.iloc[:, 0])  # Assuming a single column for the survival function
        
        SurT = list()
        SurT.append(x)
        SurT.append(y)
        SurTime.append(SurT)
        
    return SurTime

def visualRule(RiD, NBoot, NSamp, Pattern, InPutTimeStamps, Antecedent, Consequent, Confidence, Support, TID_e, SurvivalFunctions):
    dumAnt = Antecedent[RiD]
    dumCon = Consequent[RiD]

    suppAnt = list()
    dumSuppNum = 0
    for i in range(len(Pattern)):
        dum = Pattern[i]
        index = np.isin(dumAnt,dum)
        
        if np.all(index):
            suppAnt.append(i)
            dumSuppNum = dumSuppNum + 1

    suppAnt = np.array(suppAnt)


    BootConfs = list()
    BootSFunc = list()
    for i in range(NBoot):

        Candidates = np.random.permutation(len(Pattern))[:NSamp]
        Candidates = np.sort(Candidates)
        
        dumConf = Confidence[RiD]
        dumSupp = Support[RiD]
        dumTID = TID_e[RiD]
        index = np.isin(Candidates, dumTID)
        dumTID_b = Candidates[index]
        indexA = np.isin(Candidates, suppAnt)
        index = np.sum(index)
        indexA = np.sum(indexA)
        BootConfs.append(index/indexA)
        
        
        bootsfunc = calculateDistributionFunction(dumAnt,dumCon,Pattern,InPutTimeStamps,np.array(dumTID_b))
        
        bootsfuncX = SurvivalFunctions[RiD][0]
        bootsfuncY = np.interp(bootsfuncX,bootsfunc[0][0],bootsfunc[0][1])
        
        bootsfuncIn = list()
        bootsfuncIn.append(bootsfuncX)
        #bootsfuncIn.append(bootsfuncY)
        bootsfuncIn.append((1-bootsfuncY)*(index/indexA))
        
        BootSFunc.append(bootsfuncIn)
       
    SurYs = np.empty([0,len(bootsfuncIn[0])])
    for i in range(len(BootSFunc)):
        dumY = BootSFunc[i][1]
        SurYs = np.vstack([SurYs,dumY])
        
    LowerSurYs = np.percentile(SurYs, 5, axis = 0)
    UpperSurYs = np.percentile(SurYs, 95, axis = 0)

    line_width = 2

    plt.xlabel('Elapsed Time [Years]')  
    plt.ylabel('Probability')

    dum1 = str(dumAnt)
    dum2 = "⇒"
    dum3 = str(dumCon)
    title_text = f"{dum1} {dum2} {dum3}"  # Using f-string for formatting

    plt.title(title_text)

    for i in range(len(SurYs)-1):
        plt.plot(SurvivalFunctions[RiD][0], SurYs[i], color='yellow')
        
    plt.plot(SurvivalFunctions[RiD][0], SurYs[len(SurYs)-1], color='yellow', label='Bootstraps')
        
    plt.plot(SurvivalFunctions[RiD][0], (1-SurvivalFunctions[RiD][1])*Confidence[RiD], color='blue', label='Confidence function', linewidth=line_width)
    plt.plot(SurvivalFunctions[RiD][0], LowerSurYs, color='black', label='Confidence bounds', linewidth=line_width)
    plt.plot(SurvivalFunctions[RiD][0], UpperSurYs, color='black', linewidth=line_width)


    plt.legend() 
    
    return LowerSurYs, UpperSurYs

def predictEvent(TestSeq, Antecedent, Consequent, Rules, Support, Confidence):
    
    AntecedentRules = list()
    ConsequentRules = list()
    NumRules = list()
    SupportRules = list()
    ConfidenceRules = list()



    MatchedRules=np.zeros(len(Antecedent))
    for i in range(len(Antecedent)):
        dumAnt = Antecedent[i]
        indexA = np.in1d(dumAnt,TestSeq)
        
        if not np.all(indexA):
            continue
        
        dumCon = Consequent[i]
        indexA = np.in1d(dumCon,TestSeq)
        
        if np.any(indexA):
            continue
        
        MatchedRules[i] = 1 # Save RuleID
               
        dum0 = dumAnt # Save Antecedent
        AntecedentRules.append(dum0)
        
        dum01 = dumCon # Save Consequent
        ConsequentRules.append(dum01)
        
        dum2 = Rules[i] # Save Just the rule
        NumRules.append(dum2)
        
        dum4 = Support[i] # Save Support
        SupportRules.append(dum4)
        
        dum5 = Confidence[i] # Save Confidence
        ConfidenceRules.append(dum5)
        


    # After Process

    [uniqueA,uniqueB,uniqueC]=np.unique(np.array(ConsequentRules),return_index=True, return_inverse=True)
    canstay = np.empty(0)
    for i in range(len(uniqueA)):
        indexA = uniqueC == i
        foundedAnt = np.array(AntecedentRules)[indexA]
        dum = np.empty(0)
        for j in range(len(foundedAnt)):
            # len(foundedAnt[j].split("; "))
            dum = np.hstack([dum, foundedAnt[j]])
        indexB = np.max(dum,axis = 0)
        indexC = np.where(dum == indexB)[0]
        indexD = np.where(indexA)[0][indexC]

        canstay = np.hstack([canstay, indexD])

    canstay = canstay.astype(int)

    AntecedentRules = np.array(AntecedentRules)[canstay]
    ConsequentRules = np.array(ConsequentRules)[canstay]
    SupportRules = np.array(SupportRules)[canstay]
    ConfidenceRules = np.array(ConfidenceRules)[canstay]
    
    return AntecedentRules, ConsequentRules, SupportRules, ConfidenceRules