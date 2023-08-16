import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from datetime import date
from datetime import datetime
import copy

from modules.functions import *
from modules.matlab_translater import *

import joblib

# #%% Load Files


Antecedent = joblib.load('Antecedent.joblib')
Consequent = joblib.load('Consequent.joblib')
Rules = joblib.load('Rules.joblib')
Support = joblib.load('Support.joblib')
Confidence = joblib.load('Confidence.joblib')
SurvivalFunctions = joblib.load('SurvivalFunctions.joblib')
ConfidenceFunctions = joblib.load('ConfidenceFunctions.joblib')



# #%%

TestSeq = [4, 5]
TestTime = []

AntecedentRules, ConsequentRules, SupportRules, ConfidenceRules = predictEvent(TestSeq, Antecedent, Consequent, Rules, Support, Confidence)

# #%%
AfterFive = list()
for i in range(len(ConfidenceRules)):
    dumconf = ConfidenceRules[i]
    dumSx = ConfidenceFunctions[i][0]
    dumSy = ConfidenceFunctions[i][1]
    
    AfterFive.append(np.interp(5,dumSx,dumSy))
del i, dumSx, dumSy


Stack=pd.DataFrame(np.column_stack([AntecedentRules, ConsequentRules, SupportRules, ConfidenceRules, AfterFive]),columns = ["Source event", "Predicted event", "Support", "Confidence", "Probability after 5 years"])


        
