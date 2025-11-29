## 2D Coupled SMR Single Fuel Rod Analytics
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

# Constants
# ------------------------------------------------------------------------------

# Macroscopic parameters
Q = 160*1E6 # Thermal Output [MW_th -> W_th]

# 17x17 Fuel Assembly
L = 78.74*2.54 # Active Height [in -> cm]
N_Assy = 37 # Number of Fuel Assemblies
N_FP = 264 # Fuel Pins per Assy
# N_GT = 24 # Guide Tubes per Assy
# N_IT = 1 # Instrument Tubes per Assy
N_Rods = N_FP*N_Assy # Total Number of Fuel Rods
# pitch_Assy = 8.466*2.54 # Assembly Pitch [in -> cm]
pitch_Assy = 21.50 # Assembly Pitch [cm]

# Fuel Rod Data
# pitch_FP = 0.496*2.54 # Fuel Rod Pitch [in -> cm]
pitch_FP = 1.26 # Fuel Rod Pitch [cm]
# OD_FP = 0.374*2.54 # Cladding Outer Diameter [in -> cm]
# ID_FP = 0.326*2.54 # Cladding Inner Diameter [in -> cm]
c = 0.024*2.54 # Clad Thickness [in -> cm]
g = 0.0065*2.54 # Pellet-Clad Gap [in -> cm]
s = 0.3195*2.54 # Fuel Pellet Diameter [in -> cm]
R = s/2 # Fuel-pellet Radius [cm]
E = 0.0495 # Fuel Enrichment
f = 0.95 # Fuel Pellet Packing Density
rho_UO2 = 10.97 # Theoretical Density [g/cc]

# Guide/Instrument Tube Data
OD_GT = 0.482*2.54 # Guide Tube Outer Diameter [in -> cm]
ID_GT = 0.397*2.54 # Guide Tube Inner Diameter [in -> cm]

# Thermal Hydraulic Data
p = 1.2755E7 # System Pressure [Pa]
# mdot = 4.66E6*0.453592/3600 # Coolant Mass Flow Rate [lb/h->kg/s]
v = 2.7*12*2.54 # Average Coolant Velocity [ft/s -> cm/s]
T_in = (497-32)*(5/9) + 273.15 # Inlet Temperature [F->K]
T_out = (597-32)*(5/9) + 273.15 # Outlet Temperature [F->K]

# Water Properties
# mu = 9.67306535E-5/100 # Viscosity [kg/m.s] -> [kg/cm.s]
mu = 0.01 # Viscosity [kg/cm.s]
C_p = 5066.028 # Isobaric Heat Capacity [J/kg.K]
k_H2O = 0.005936 # Thermal Conductivity [W/cm.K]
rho_H2O = 0.765202/1000 # Density [g/cc -> kg/cc]

# AP1000 Supplemental Properties
P_ci = 1000 # contact pressure between ci and fuel pellet surface [psia]
k_c = 10.5*0.017295772056 # Cladding Thermal Conductivity [BTU/h.ft.R] -> [W/cm.K]
alpha = 0.6 # gap heat transfer pressure coefficient

# Calculated Properties
V_fuel = N_Rods*np.pi*L*R**2 # Total Core Fuel Volume [cc]
q3prime = Q/V_fuel; # Average Fuel Volumetric Heat Rate [W/cc]
q_s = Q/N_Rods # Total Heat Generation in a Single Fuel Element [W]
# T_infty = (T_in+T_out)/2 # Bulk Fluid Temperature [K]
T_infty = (543-32)*(5/9) + 273.15 # Avg Coolant Temperature [F->K]

# CASMO Input Params 
PDE = Q/(L*np.pi*np.power(R,2)*N_Rods) # [kW/L]
print('PDE', PDE)
PRE = p*1e-5 # [bar]
print('PRE', PRE)
print('PIN 1  ', R, ' ', R+g, ' ', R+g+c, '* FUEL PIN')
print('PIN 2 ', ID_GT/2, ' ', OD_GT/2, ' /  \'COO\'    \'BOX\'  * INSTR TUBE') 
print('PIN 3 ', ID_GT/2, ' ', OD_GT/2, '/   \'COO\'    \'BOX\'  * GUIDE TUBE')
print('FUE 1,', rho_UO2*f, '/', E*1E2)
BORs = [0.0, 600.0, 1200.0]
TMs = [531.5, 557.0, 587.0]
TFs = [600.0, 850.0, 1200.0]
i=0
for BOR in BORs:
    for TM in TMs:
        for TF in TFs:
            i+=1
            print('TTL *+ PERTURBATION', i)
            print('TFU', TF)
            print('BOR', BOR)
            print('TMO', TM)
            print('AVE \'CELL\', / 0 1 0 0 0 0')
print()
i=0
for BOR in BORs:
    for TM in TMs:
        for TF in TFs:
            i+=1
            print('TTL *+ PERTURBATION', i)
            print('TFU', TF)
            print('BOR', BOR)
            print('TMO', TM)
            print('REF',pitch_Assy,'/\'MOD\' / / /')
            print('AVE \'CELL\', / 0 1 0 0 0 0')

# Computation Parameters
Nr = 50 # discrete radial cells
convergence_tolerance = 0.001 # convergence tolerance
# sublayerSpace = np.linspace(R+g+c,R_FC,Nr)
claddingSpace = np.linspace(R+g,R+g+c,Nr)
gapSpace = np.linspace(R,R+g,Nr)
fuelSpace = np.linspace(0,R,Nr)

# Load Thermal Data
# ------------------------------------------------------------------------------

class InterpolatedProperties:
    def __init__(self, rhoC_pCSV, kCSV):
        self.T1 = []
        self.rho = []
        self.cp = []
        self.T2 = []
        self.kt = []

        with open(rhoC_pCSV, 'r') as f:
            reader = csv.reader(f)
            next(reader)  # skip header
            for row in reader:
                T_val = float(row[0])
                rho_val = float(row[1])
                cp_val = float(row[2])

                self.T1.append(T_val)
                self.rho.append(rho_val)
                self.cp.append(cp_val)

        with open(kCSV, 'r') as f:
            reader = csv.reader(f)
            next(reader)  # skip header
            for row in reader:
                T_val = float(row[0])
                k_val = float(row[1])

                self.T2.append(T_val)
                self.kt.append(k_val)

        # Convert to arrays
        self.T1 = np.array(self.T1)
        self.rho = np.array(self.rho)
        self.cp = np.array(self.cp)
        self.T2 = np.array(self.T2)
        self.kt = np.array(self.kt)

        # Create interpolators
        self.rho_interp = interp1d(self.T1, self.rho, kind='linear', fill_value='extrapolate')
        self.cp_interp = interp1d(self.T1, self.cp, kind='linear', fill_value='extrapolate')
        self.k_interp = interp1d(self.T2, self.kt, kind='linear', fill_value='extrapolate')

    def rhoc_p(self, T, f):
        rho_val = float(self.rho_interp(T))
        cp_val = float(self.cp_interp(T))
        return f*rho_val*cp_val
        
    def k(self, T):
        k_val = float(self.k_interp(T))
        return k_val


UO2_rhoC_p = 'proj6Data\\thermalData\\UO2_rhoC_p.csv'
UO2_k = 'proj6Data\\thermalData\\UO2_k.csv'
fuelObj = InterpolatedProperties(UO2_rhoC_p, UO2_k)

# Functions
# ------------------------------------------------------------------------------

# Cladding Temperature Profile
# Accepts total heat generation q_s, radius r (between R+g+c and R+g)
#         and cladding outside temperature T_co
# Returns temperature at radius r
def Tclad(q_s, r, T_co):
    T = T_co + (q_s*np.log((R+g+c)/r))/(2*np.pi*L*k_c)
    return T

# Gap Temperature Profile
# Gap heat transfer is dependent on temperature, so substitutes temperatures 
# back in to converge on true value
# Accepts total heat generation q_s, radius r (between R+g and R)
#         and cladding inside temperature T_ci
# Returns temperature at radius r
def Tgap(q_s, r, T_ci):
    T = 1400
    T_guess = 0
    while np.abs(T_guess-T)/T > convergence_tolerance:
        # print(T) # debug
        T_guess = T
        h_gap = hGap(T_guess, T_ci)
        T = T_ci + (q_s*(1+((np.log(r/R))/(np.log(R/(R+g))))))/(2*np.pi*L*R*h_gap)
    return T

# Fuel Temperature Profile
# Fuel heat transfer is dependent on temperature, so substitutes temperatures 
# back in to converge on true value
# Accepts total heat generation q_s, radius r (between R and 0)
#         and fuel surface temperature T_s
# Returns temperature at radius r
def Tfuel(q_s, r, T_s):
    T = 2000
    T_guess = 0
    while np.abs(T_guess-T)/T > convergence_tolerance:
        # print(T) # debug
        T_guess = T
        k_f = fuelObj.k(T_guess)
        T = T_s + (q_s*(R**2-r**2))/(L*np.pi*(R**2)*4*k_f)
    return T

# gap Heat Transfer
# Accepts Fuel pellet surface temp (T_s) and Cladding inside surface temp
#         (T_ci) in Rankine
# Returns gap heat transfer coefficient h_gap (BTU/h-ft^2-R) -> [W/cm^2.K]
def hGap(T_s, T_ci):
    T_gas = ((T_s+T_ci)/2)*9/5 # [K->R]
    R_ci = (R+g)/30.48 # [cm->ft]
    R_s = R/30.48
    h_gap = alpha*P_ci + kGas(T_gas)/((R_ci-R_s)+14.4e-6)
    return h_gap*5.678263398/10000
    
# Helium Conductivity
# Accepts gas temperature T_gas in Rankine
# Returns gas thermal conductivity k_gas (BTU/h-ft-R)
def kGas(T_gas):
    k_gas = 1.314e-4*(T_gas**0.668)
    return k_gas

# Fuel Conductivity
# Accepts fuel temperature T_f in Rankine
# Returns fuel conductivity k_f (BTU/h-ft-R)
def kFuel(T_f):
    T = T_f
    k = (1/(11.8 + (0.0238*T))) + ((8.775e-13)*(T**3))
    return k

# Script
# ------------------------------------------------------------------------------

# Bulk Fluid Temperature Value
# -----------------------------------------------------------------------------

print()
print("T_infty =", str.format('{0:.2f}', T_infty), "K") # print
print()

# Laminar Sublayer Temperature Distribution
# -----------------------------------------------------------------------------

h_coolant = 35000/10000 # Water Convective Heat Transfer coeff. [W/m^2.K -> W/cm^2.K]

h_coolant = rho_H2O*C_p*v

T_co = T_infty + q_s/(2*np.pi*(R+g+c)*L*h_coolant) # Cladding Outer Temp [K]
    
print()
print("T_co =", str.format('{0:.2f}', T_co), "K") # print
print()

# Cladding Temperature Distribution
# -----------------------------------------------------------------------------
    
Tcladr = np.zeros(Nr) # init
for i in range(Nr):
    Tcladr[i] = Tclad(q_s, claddingSpace[i], T_co)
        
T_ci = np.max(Tcladr) # cladding outside temp
print()
print("T_ci =", str.format('{0:.2f}', T_ci), "K") # print
print()
    
plt.plot(claddingSpace, Tcladr) # plot radial temp profile

# Gas Gap Temperature Distribution
# -----------------------------------------------------------------------------
    
Tgapr = np.zeros(Nr) # init
for i in range(Nr):
    Tgapr[i] = Tgap(q_s, gapSpace[i], T_ci)
        
T_s = np.max(Tgapr) # cladding outside temp
print()
print("T_s =", str.format('{0:.2f}', T_s), "K") # print
print()
    
plt.plot(gapSpace, Tgapr) # plot radial temp profile

# Fuel Temperature Distribution
# -----------------------------------------------------------------------------
    
Tfuelr = np.zeros(Nr) # init
for i in range(Nr):
    Tfuelr[i] = Tfuel(q_s, fuelSpace[i], T_s)
        
T_m = np.max(Tfuelr) # cladding outside temp
print()
print("T_c =", str.format('{0:.2f}', T_m), "K") # print
print()
    
plt.plot(fuelSpace, Tfuelr) # plot radial temp profile

# Finalize Plot
# -----------------------------------------------------------------------------
plt.axvline(x = 0, color = 'k', label = 'Fuel Centerline', alpha=0.3, linewidth=1)
plt.axvline(x = R, color = 'k', label = 'Fuel Surface', alpha=0.3, linewidth=1)
plt.axvline(x = R+g, color = 'k', label = 'Cladding Inside Surface', alpha=0.3, linewidth=1)
plt.axvline(x = R+g+c, color = 'k', label = 'Cladding Outside Surface', alpha=0.3, linewidth=1)
plt.xlabel('Radius (cm)')
plt.ylabel('Temperature (K)')
plt.title('Radial Temperature Profile')
plt.tight_layout()
plt.savefig('results\\project6\\radialTemperatureCurve.jpg')
plt.show()
plt.close()

k_favg = np.average(kFuel(Tfuelr)) # RADIALLY averaged fuel conductivity [W/cm.K]
T_favg = np.average(Tfuelr) # Radial average fuel temperature [K]
rhoC_pavg = fuelObj.rhoc_p(T_favg, f) # Approximate volumetric heat capacity [J/cm^3.K]
h_gapavg = np.average(hGap(Tgapr, T_ci)) # Radial average gap heat transfer [W/cm^2.K]
k_gapavg = -h_gapavg*np.log(R/(R+g))*(R+g/2) # Equivalent gap conductivity [W/cm.K]

print("k_f =", str.format('{0:.5f}', k_favg), "W/cm.K") # print
print("rhoC_p_f =", str.format('{0:.5f}', rhoC_pavg), "J/cm^3.K")
print("k_gap =", str.format('{0:.5f}', k_gapavg), "W/cm.K")
print("k_c =", str.format('{0:.5f}', k_c), "W/cm.K")
print("rhoC_p_coolant =", str.format('{0:.5f}', rho_H2O*C_p), "J/cm^3.K") # print
print("h_coolant =", str.format('{0:.5f}', h_coolant), "W/cm^2.K") # print
print("q3prime =", str.format('{0:.5f}', q3prime), "W/cm^3")
print("q_s =", str.format('{0:.5f}', q_s), "W")
print()

# CSV-ify
# k_FA = 8*np.pi*k_favg
k_FA = 1/( (1/(8*np.pi*k_favg)) + (np.log((R+g)/R)/(2*np.pi*k_gapavg)) + (np.log((R+g+c)/(R+g))/(2*np.pi*k_c)) + (1/(2*np.pi*(R+g+c)*h_coolant)) ) # Effective Thermal Conductivity of a Fuel Assembly [W/cm.K]

FA_CSV = 'proj6Data\\thermalData\\Fuel_Lump.csv'
header = ["k", "rhoCp", "conv"]
row = [k_FA, rhoC_pavg, 0]
with open(FA_CSV, "w", newline="") as myFile:
    writer = csv.writer(myFile)
    writer.writerow(header)
    writer.writerow(row)

MOD_CSV = 'proj6Data\\thermalData\\H2O_Lump.csv'
header = ["k", "rhoCp", "conv"]
row = [0.006, rho_H2O*C_p, 1]
with open(MOD_CSV, "w", newline="") as myFile:
    writer = csv.writer(myFile)
    writer.writerow(header)
    writer.writerow(row)


# -----------------------------------------------------------------------------
# Redo with Collapsed Properties
# -----------------------------------------------------------------------------

# Cladding Temperature Profile
# Accepts total heat generation q_s, radius r (between R+g+c and R+g)
#         and cladding outside temperature T_co
# Returns temperature at radius r
def Tclad2(q_s, r, T_co, k_c):
    T = T_co + (q_s*np.log((R+g+c)/r))/(2*np.pi*L*k_c)
    return T

# Gap Temperature Profile
# Accepts total heat generation q_s, radius r (between R+g and R)
#         and cladding inside temperature T_ci
# Returns temperature at radius r
def Tgap2(q_s, r, T_ci, k_gap):
    T = T_ci + (q_s*np.log((R+g)/r))/(2*np.pi*L*k_gap)
    return T

# Fuel Temperature Profile
# Fuel heat transfer is dependent on temperature, so substitutes temperatures 
# back in to converge on true value
# Accepts total heat generation q_s, radius r (between R and 0)
#         and fuel surface temperature T_s
# Returns temperature at radius r
def Tfuel2(q_s, r, T_s, k_f):
    T = T_s + (q_s*(R**2-r**2))/(L*np.pi*(R**2)*4*k_f)
    return T

def fullDomain(r, q_s, T_co, T_ci, T_s, k_c, k_gap, k_f):
    if r < R:
        return Tfuel2(q_s, r, T_s, k_f)
    elif r < R + g:
        return Tgap2(q_s, r, T_ci, k_gap)
    elif r < R + g + c:
        return Tclad(q_s, r, T_co)
    return T_infty

# Bulk Fluid Temperature Value
# -----------------------------------------------------------------------------

print()
print("T_infty =", str.format('{0:.2f}', T_infty), "K") # print
print()

# Laminar Sublayer Temperature Distribution
# -----------------------------------------------------------------------------

# D_eq = 4*(pitch**2 - np.pi*(myOD/2)**2)/(2*np.pi*(myOD/2)) # Hydraulic Diameter [cm]
# Re = D_eq*mdot/mu # Reynolds Number
# Pr = C_p*mu/k_H2O # Prandtl Number
# h_coolant = (0.023*k_H2O/D_eq)*Re**0.8*Pr**0.4 # Water Convective Heat Transfer Coefficient via Dittus-Boelter [W/m^2.K]
h_coolant = 3.5

T_co = T_infty + q_s/(2*np.pi*(R+g+c)*L*h_coolant) # Cladding Outer Temp [K]
    
print()
print("T_co =", str.format('{0:.2f}', T_co), "K") # print
print()

# Cladding Temperature Distribution
# -----------------------------------------------------------------------------
    
Tcladr = np.zeros(Nr) # init
for i in range(Nr):
    Tcladr[i] = Tclad2(q_s, claddingSpace[i], T_co, k_c)
        
T_ci = np.max(Tcladr) # cladding outside temp
print()
print("T_ci =", str.format('{0:.2f}', T_ci), "K") # print
print()
    
plt.plot(claddingSpace, Tcladr) # plot radial temp profile

# Gas Gap Temperature Distribution
# -----------------------------------------------------------------------------
    
Tgapr = np.zeros(Nr) # init
for i in range(Nr):
    Tgapr[i] = Tgap2(q_s, gapSpace[i], T_ci, k_gapavg)
        
T_s = np.max(Tgapr) # cladding outside temp
print()
print("T_s =", str.format('{0:.2f}', T_s), "K") # print
print()
    
plt.plot(gapSpace, Tgapr) # plot radial temp profile

# Fuel Temperature Distribution
# -----------------------------------------------------------------------------
    
Tfuelr = np.zeros(Nr) # init
for i in range(Nr):
    Tfuelr[i] = Tfuel2(q_s, fuelSpace[i], T_s, k_favg)
        
T_m = np.max(Tfuelr) # cladding outside temp
print()
print("T_c =", str.format('{0:.2f}', T_m), "K") # print
print()
    
plt.plot(fuelSpace, Tfuelr) # plot radial temp profile

# Finalize Plot
# -----------------------------------------------------------------------------
plt.axvline(x = 0, color = 'k', label = 'Fuel Centerline', alpha=0.3, linewidth=1)
plt.axvline(x = R, color = 'k', label = 'Fuel Surface', alpha=0.3, linewidth=1)
plt.axvline(x = R+g, color = 'k', label = 'Cladding Inside Surface', alpha=0.3, linewidth=1)
plt.axvline(x = R+g+c, color = 'k', label = 'Cladding Outside Surface', alpha=0.3, linewidth=1)
plt.xlabel('Radius (cm)')
plt.ylabel('Temperature (K)')
plt.title('Radial Temperature Profile')
plt.tight_layout()
plt.savefig('results\\project6\\radialTemperatureCurve2.jpg')
plt.show()
plt.close()
