import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

from modules.functions import *


import joblib
InPutSequences = joblib.load('InPutSequences.joblib') # a -> b <==> [a -1 b]
InPutTimeStamps = joblib.load('InPutTimeStamps.joblib')
    

# #%% Writing input file for sequential rule miner
Pattern=InPutSequences

writeSPMF(Pattern)

# #%% Exectuting sequential rule mining 

import os

command =["java -jar spmf.jar run RuleGrowth input.txt output.txt 10% 10%"][0]
#command =["java -jar spmf.jar run TopSeqRules input.txt output.txt 100 10%"][0]
os.system(command)

# #%% Read frequent sequences


Rules, Support, Confidence = readSeqRules()

# #%% Finding supporting Trace ID-s of rules

TID_e = findTID(Rules, Pattern)



# #%% Making Antecedent and Consequent Vectors

Antecedent, Consequent = separateRule(Rules)

# #%% Calculated the distribution of elapsed times with Kaplan Meier estimator

SurvivalFunctions = calculateDistributionFunctions(Antecedent,Consequent,Pattern,InPutTimeStamps,TID_e)

# #%% Calculate the time dependent confidence functions

ConfidenceFunctions = calculateConfidenceFunctions(SurvivalFunctions,Confidence)

# #%% Estimating Confidence Intervals by bootstrapping and plot results for one rule

RiD = 30 # ID of rule
NBoot = 200 # Number of Bootstraps
NSamp = 200 # Number of Samples


LowerSurYs, UpperSurYs = visualRule(RiD, NBoot, NSamp, Pattern, InPutTimeStamps, Antecedent, Consequent, Confidence, Support, TID_e, SurvivalFunctions)


# #%% Write supplementary Variables

joblib.dump(Antecedent, 'Antecedent.joblib')
joblib.dump(Consequent, 'Consequent.joblib')
joblib.dump(Rules, 'Rules.joblib')
joblib.dump(Support, 'Support.joblib')
joblib.dump(Confidence, 'Confidence.joblib')
joblib.dump(SurvivalFunctions, 'SurvivalFunctions.joblib')
joblib.dump(ConfidenceFunctions, 'ConfidenceFunctions.joblib')





    




