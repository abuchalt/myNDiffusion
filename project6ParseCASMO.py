## CASMO Parsing Script
# ------------------------------------------------------------------------------
# Adam Buchalter
#
# desc
# ------------------------------------------------------------------------------

# Packages
# ------------------------------------------------------------------------------
import numpy as np
import matplotlib.pyplot as plt
import scipy.io as sio
from scipy.interpolate import interp1d
import os
import csv
import re
import pandas as pd

# Functions
# ------------------------------------------------------------------------------
def parseCASMOOutput(filepath, csvout):

    dictList = []

    with open(filepath, 'r') as f:
        lines = f.readlines()
    
    for i in range(len(lines)):

        currentLine = lines[i]

        if "CELL AVERAGED CROSS SECTIONS" in currentLine:

            myDict = {}

            match = re.findall(r"[-+]?\d*\.?\d+", lines[i-4])
            myDict['BURNUP'] = float(match[0])
            # myDict['V'] = float(match[1])
            myDict['TF'] = float(match[2])
            myDict['TM'] = float(match[3])
            myDict['BOR'] = float(match[4])

            myDict['group'] = []
            myDict['nuSigma_f'] = []
            myDict['Sigma_f'] = []
            myDict['Sigma_a'] = []
            myDict['D'] = []
            myDict['Sigma_R'] = []
            myDict['Chi'] = []
            myDict['Sigma_s1g'] = []
            myDict['Sigma_s2g'] = []

            match = re.findall(r"[-+]?\d+\.?\d*[eE][-+]?\d+", lines[i+4])
            group = 1
            nuSigma_f = float(match[3])
            Sigma_f = float(match[2])
            Sigma_a = float(match[0])
            D = float(match[6])
            Sigma_R = Sigma_a + float(match[5])
            Chi = 1
            match = re.findall(r"[-+]?\d+\.?\d*[eE][-+]?\d+", lines[i+9])
            Sigma_s1g = float(match[0])
            match = re.findall(r"[-+]?\d+\.?\d*[eE][-+]?\d+", lines[i+11])
            Sigma_s2g = float(match[0])
            myDict['group'].append(group)
            myDict['nuSigma_f'].append(nuSigma_f)
            myDict['Sigma_f'].append(Sigma_f)
            myDict['Sigma_a'].append(Sigma_a)
            myDict['D'].append(D)
            myDict['Sigma_R'].append(Sigma_R)
            myDict['Chi'].append(Chi)
            myDict['Sigma_s1g'].append(Sigma_s1g)
            myDict['Sigma_s2g'].append(Sigma_s2g)

            match = re.findall(r"[-+]?\d+\.?\d*[eE][-+]?\d+", lines[i+5])
            group = 2
            nuSigma_f = float(match[3])
            Sigma_f = float(match[2])
            Sigma_a = float(match[0])
            D = float(match[6])
            Sigma_R = Sigma_a + float(match[5])
            Chi = 0
            match = re.findall(r"[-+]?\d+\.?\d*[eE][-+]?\d+", lines[i+9])
            Sigma_s1g = float(match[1])
            match = re.findall(r"[-+]?\d+\.?\d*[eE][-+]?\d+", lines[i+11])
            Sigma_s2g = float(match[1])
            myDict['group'].append(group)
            myDict['nuSigma_f'].append(nuSigma_f)
            myDict['Sigma_f'].append(Sigma_f)
            myDict['Sigma_a'].append(Sigma_a)
            myDict['D'].append(D)
            myDict['Sigma_R'].append(Sigma_R)
            myDict['Chi'].append(Chi)
            myDict['Sigma_s1g'].append(Sigma_s1g)
            myDict['Sigma_s2g'].append(Sigma_s2g)

            dictList.append(myDict)
        
    # Determine the fieldnames (CSV headers)
    # This ensures all possible keys from all dictionaries are included
    # and handles cases where dictionaries might have different keys.
    fieldnames = sorted(list(set(key for d in dictList for key in d.keys())))

    # Open the CSV file in write mode
    with open(csvout, 'w', newline='') as csvfile:
        # Create a DictWriter object
        writer = csv.DictWriter(csvfile, fieldnames=fieldnames)

        # Write the header row
        writer.writeheader()

        # Write all the dictionaries as rows
        writer.writerows(dictList)

    print(f"Data successfully written to {csvout}")
    return dictList

myCASMO = 'proj6Data\\CASMO\\nuscalePwrREF.out'
myCSV = 'proj6Data\\neutronData\\nuscalePwrREF.csv'
myMat = 'proj6Data\\neutronData\\interpTablesREF.mat'

myDictList = parseCASMOOutput(myCASMO, myCSV)

# BORs = [0.0, 200.0, 400.0, 600.0, 800.0, 1000.0]
BORs = [1000.0, 2500.0, 4000.0]
# BURNUPs = [0.0, 0.1, 0.2, 0.5, 1.0, 2.0, 5.0, 10.0, 20.0, 50.0]
# BURNUPS = [0.0]
# TMs = [532.04, 547.04, 552.04, 557.04, 562.04, 567.04, 582.04]
TMs = [500, 557.0, 620.0]
# TFs = [650.0, 750.0, 850.0, 950.0, 1050.0, 1150.0]
TFs = [600.0, 850.0, 1200.0]

Chi1 = 1
Chi2 = 0
D1 = np.zeros((len(TFs), len(TMs), len(BORs)))
D2 = np.zeros_like(D1)
Sigma_R1 = np.zeros_like(D1)
Sigma_R2 = np.zeros_like(D1)
Sigma_a1 = np.zeros_like(D1)
Sigma_a2 = np.zeros_like(D1)
# Sigma_f = np.zeros_like(D1)
# Sigma_f = np.zeros_like(D1)
Sigma_s11 = np.zeros_like(D1)
Sigma_s12 = np.zeros_like(D1)
Sigma_s21 = np.zeros_like(D1)
Sigma_s22 = np.zeros_like(D1)
nuSigma_f1 = np.zeros_like(D1)
nuSigma_f2 = np.zeros_like(D1)

g1ind = 0
g2ind = 1

# Convert to Format MATLAB likes for interpolation
for myDict in myDictList:
    
    TF_ind = TFs.index(myDict['TF'])
    TM_ind = TMs.index(myDict['TM'])
    BOR_ind = BORs.index(myDict['BOR'])
    # BURNUP_ind = BURNUPs.index(myDict['BURNUP'])
    
    D1[TF_ind, TM_ind, BOR_ind] = myDict['D'][g1ind]
    D2[TF_ind, TM_ind, BOR_ind] = myDict['D'][g2ind]
    Sigma_R1[TF_ind, TM_ind, BOR_ind] = myDict['Sigma_R'][g1ind]
    Sigma_R2[TF_ind, TM_ind, BOR_ind] = myDict['Sigma_R'][g2ind]
    Sigma_a1[TF_ind, TM_ind, BOR_ind] = myDict['Sigma_a'][g1ind]
    Sigma_a2[TF_ind, TM_ind, BOR_ind] = myDict['Sigma_a'][g2ind]
    Sigma_s11[TF_ind, TM_ind, BOR_ind] = myDict['Sigma_s1g'][g1ind]
    Sigma_s12[TF_ind, TM_ind, BOR_ind] = myDict['Sigma_s1g'][g2ind]
    Sigma_s21[TF_ind, TM_ind, BOR_ind] = myDict['Sigma_s2g'][g1ind]
    Sigma_s22[TF_ind, TM_ind, BOR_ind] = myDict['Sigma_s2g'][g2ind]
    nuSigma_f1[TF_ind, TM_ind, BOR_ind] = myDict['nuSigma_f'][g1ind]
    nuSigma_f2[TF_ind, TM_ind, BOR_ind] = myDict['nuSigma_f'][g2ind]

data_to_save = {
    'TFs': TFs,
    'TMs': TMs,
    'BORs': BORs,
    'D1': D1,
    'D2': D2,
    'Sigma_R1': Sigma_R1,
    'Sigma_R2': Sigma_R2,
    'Sigma_a1': Sigma_a1,
    'Sigma_a2': Sigma_a2,
    'Sigma_s11': Sigma_s11,
    'Sigma_s12': Sigma_s12,
    'Sigma_s21': Sigma_s21,
    'Sigma_s22': Sigma_s22,
    'nuSigma_f1': nuSigma_f1,
    'nuSigma_f2': nuSigma_f2
}

sio.savemat(myMat, data_to_save)