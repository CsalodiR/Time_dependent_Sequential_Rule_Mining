import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

from modules.functions import *
import joblib


#Data is synthetised with Markov models

NoEvents = 10 # No of unique events
SeqEndProb = 0.4 # Probability of ending the sequence
StampEndProb = 0.5 # Probability of event is handled together
NoInPuts = 500 # No of samples



StateTransition = np.random.rand(NoEvents,NoEvents)
StateTransition = StateTransition / np.sum(StateTransition,axis = 0) # Making state transition matrix


InPutSequences = list()



for i in range(NoInPuts): # Generate events with roulette wheel selection
    
    actseq = np.array([])
    firstEvent = np.floor(np.random.rand(1)*NoEvents)
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
            
            actseq = np.concatenate([actseq, np.array([-2])])
            makeMore = False
        else:
            actseq = np.concatenate([actseq, np.array([-1])])
        j=j+2
    actseq = actseq.astype(int)
    InPutSequences.append(actseq)

#%%
# Parameters of event time distribution (exponential)
# Timestamps are provided in datenum format
# Initial date of sequence is between 2005 (732313) and 2008 (733408)

initL = 732313
initU = 733408

# Parameters of exponential distribution of elapsed times
ExponentialParameters = np.random.randint(50, 100 + 1,size=NoEvents)

startDate = np.random.randint(initL, initU + 1)
np.random.exponential(scale=100.0, size=None)

InPutTimeStamps = list()
for i in range(NoInPuts):
    dumSeq = InPutSequences[i]
    times = list()
    
    startDate = np.random.randint(initL, initU + 1)
    times.append(startDate)
    actDate = startDate
    for j in range(1,len(dumSeq)):
        if dumSeq[j] == -1:
            times.append(int(-1))
            actDate = actDate + np.ceil(np.random.exponential(scale=ExponentialParameters[int(dumSeq[j-1])], size=None))
        else:
            times.append(actDate)
    times.append(int(-2))
    times = np.array(times)
    times = times.astype(int)
    InPutTimeStamps.append(times)

# Save variables
joblib.dump(InPutSequences, 'InPutSequences.joblib')
joblib.dump(InPutTimeStamps, 'InPutTimeStamps.joblib')


