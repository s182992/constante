import sympy
import math
import numpy
import matplotlib.pyplot as plt
sympy.__version__
from sympy import *
from math import log10


pK_h2co3_stdrd = 6.351 # valeur de pK a 25°C (temperature standard)
exch_H_h2co3 = 9.15    # en KJ/mol
exch_Cp_h2co3 = -371.0   # en J/(K.mol)
pK_hco3_stdrd = 10.329
exch_H_hco3 = 14.7  # en KJ/mol
exch_Cp_hco3 = -249.0  # en J/(K.mol)
pK_h3L_stdrd = 3.128
exch_H_h3L = 4.07  # en KJ/mol
exch_Cp_h3L = -131.0  # en J/(K.mol)
pK_h2L_stdrd = 4.761
exch_H_h2L = 2.28  # en KJ/mol
exch_Cp_h2L = -178.0  # en J/(K.mol)
pK_hL_stdrd = 6.396
exch_H_hL = -3.38  # en KJ/mol
exch_Cp_hL = -254.0  # en J/(K.mol)
pKe_stdrd = 14.0  # produit ionique de l'eau
a = 6894.0         # constante dans l'expression du produit ionique
pK_CaCO3_stdrd = 10**(6.2**(-9))
exch_H_CaCO3 = -1207.4  # en KJ/mol
exch_Cp_CaCO3 = 81.9  # en J/(K.mol)
exch_H_R = 2400 # en Kelvin
#H_stdrd = 3.3*10**(-4) # en mol/(m^3.Pa)
H_stdrd = 0.033 # en mol/(L.atm)
T_stdrd = 298   # température standard en Kelvin
R = 8.314472     # constante de gaz parfait en J/(K.mol)
ln_10 = 2.303    # ln(a)= ln(10)*log(a)

# constantes dans l'expression de la constante de solubilité/précipitation
A = 7.8156
B = 0.03111
C = 1502
D = 5.518

# Déclaration des variables
K_h2co3, K_hco3, K_AH3AH2, K_AH2AH, K_AHA, Kw, K_CaCO3, Hcp_CO2 = symbols('K_h2co3 K_hco3 K_AH3AH2 K_AH2AH K_AHA Kw K_CaCO3 Hcp_CO2')
T = input("T= ")  # permet d'entrer la température

# David R. Lide - CRC Handbook of Chemistry and Physics, 84th Edition-CRC Press (2003). page 1193
equations = [
    K_h2co3 - 10**(-((1/(R * ln_10)) * (((ln_10 * R * int(T) * pK_h2co3_stdrd)/T_stdrd) - exch_H_h2co3 * ((1/T_stdrd) - (1/int(T))) - exch_Cp_h2co3 * ((T_stdrd/int(T)) - 1 + ln_10 * log10(int(T)/T_stdrd))))),
    K_hco3 - 10**(-(1/(R * ln_10)) * (((ln_10 * R * int(T) * pK_hco3_stdrd)/T_stdrd) - exch_H_hco3 * ((1/T_stdrd) - (1/int(T))) - exch_Cp_hco3 * ((T_stdrd/int(T)) - 1 + ln_10 * log10(int(T)/T_stdrd)))),
    K_AH3AH2 - 10**(-(1/(R * ln_10)) * (((ln_10 * R * int(T) * pK_h3L_stdrd)/T_stdrd) - exch_H_h3L * ((1/T_stdrd) - (1/int(T))) - exch_Cp_h3L * ((T_stdrd/int(T)) - 1 + ln_10 * log10(int(T)/T_stdrd)))),
    K_AH2AH - 10**(-(1/(R * ln_10)) * (((ln_10 * R * int(T) * pK_h2L_stdrd)/T_stdrd) - exch_H_h2L * ((1/T_stdrd)-(1/int(T))) - exch_Cp_h2L * ((T_stdrd/int(T)) - 1 + ln_10 * log10(int(T)/T_stdrd)))),
    K_AHA - 10**(-(1/(R * ln_10)) * (((ln_10 * R * int(T) * pK_hL_stdrd)/T_stdrd) - exch_H_hL * ((1/T_stdrd)-(1/int(T))) - exch_Cp_hL * ((T_stdrd/int(T)) - 1 + ln_10 * log10(int(T)/T_stdrd)))),

    Kw - 10**(-(((a/ln_10)*((1/int(T))-(1/T_stdrd))) + pKe_stdrd)),

    #K_CaCO3 - 10**(-(1/(R * ln_10)) * (((ln_10 * R * int(T) * pK_CaCO3_stdrd)/T_stdrd) - exch_H_CaCO3 * ((1/T_stdrd)-(1/int(T))) - exch_Cp_CaCO3 * ((T_stdrd/int(T)) - 1 + ln_10 * log10(int(T)/T_stdrd)))),
    K_CaCO3 - 10**(-(A + B * int(T) + C/int(T) - D * log10(int(T)))),

    Hcp_CO2 - 10**(((exch_H_R/ln_10) * ((-1/int(T)) + (1/T_stdrd))) + log10(H_stdrd)),
]
A = solve(equations, K_h2co3, K_hco3, K_AH3AH2, K_AH2AH, K_AHA, Kw, K_CaCO3, Hcp_CO2)

#print(pK_h2co3, pK_hco3, pK_h3L, pK_h2L, pK_hL, pKe)

print(A)
