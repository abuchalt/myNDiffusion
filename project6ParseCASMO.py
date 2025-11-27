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
            myDict['V'] = float(match[1])
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

myCASMO = 'proj6Data\\CASMO\\pwrExample.out'
myCSV = 'proj6Data\\neutronData\\test.csv'

parseCASMOOutput(myCASMO, myCSV)