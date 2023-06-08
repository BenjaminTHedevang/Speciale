# -*- coding: utf-8 -*-
"""
Created on Wed Oct  5 08:45:12 2022
@author: Benjamin Terndrup Hedevang & Frederik Bang Sørensen
"""

# Modul import
import numpy as np

Parameter = 'Poulsen' # Navn på ændret parameter
Default = ''
#Model:
dt = 1/480 # Tidsskridt, [min]
dz = 1/4 # vertikal cellestørrelse [cm]
regn_varighed = 20 # min
regn_intensitet = 0.7 #cm/min
Poulsen = True # Brug af Ksat regnet ved Poulsen (True) eller målt (False)
Tylle = False # Benyttelse af jordegenskaber fra Tylstrup (52 - 76 cm)
g_dir = 1 # Retning af gravitation, 1 = Ned, 0 = Horisontal, -1 = Op (Kapillær)
rhob = 1.7 # Øget kompaktering [g/cm3]
ror_celle = 60 # Celle nummer hvorfra øget kompaktering sker
K_s = 19.61/60 # cm/min
#Jordegenskaber (skrueparametre)
b = 2.28  # Campbells, b [-], Gennemsnit af Campbell fra vG fit 
psi_entry = -2.19  # Air-entry pressure [cm H2O], Gennemsnit af Campbell fra vG fit (Tabel 3.3, Fig 3.5)
theta_s = 0.434  # Water content at saturation [cm^3 H2O / cm^3 soil] (pF0)
theta_wp = 0.016694  # Water content at wilting point [cm^3 / cm^3] (pF4.2)
theta_airdry = 0.002827


GW_dist = 35 # afstand til grundvandsspejl, her i sandboksen, [cm]
stepsave = 2  # Gem hver nth step
n_cells = int(GW_dist/dz)  # Antal celler
depth = np.array(np.arange(dz/2, dz*n_cells, dz))  # Midten af modelcellerne (Dybde) [cm]
t_n = int(dt**(-1)*60*2)   # Antallet af tidsskridt [-]



# Nedbør
regn_array = np.zeros(t_n)
if g_dir == 1:
    regn_array[:int(dt**(-1)*regn_varighed)] = regn_intensitet

#Initialbetingelser:
theta_hot = np.load('Hot_start72t025cm025s.npy')[-1,:] if Tylle == False else np.load('Jordfugtighed_hotstartTylle.npy')[-1,:]
F_acc = [0] # Akkumuleret infiltration
rhoref = 1.5 # Normal kompaktering
tt = [regn_array[0]*dt/dz] #Top reservoir liste
tr = regn_array[0]*dt/dz # Initial condition of top reservoir

# Infiltrationsparametre: 
GI_idle = np.array([0.11, 0.13, 0.14, 0.16, 0.17, 0.19, 0.20, 0.22,
                    0.23, 0.25, 0.27, 0.30, 0.33, 0.36, 0.40, 0.44,
                    0.52, 0.63, 0.76, 0.87, 0.96, 1.00, 1.00, 1.00,
                    1.00, 1.00, 1.00, 0.97, 0.94, 0.87, 0.84, 0.78,
                    0.72, 0.67, 0.61, 0.56, 0.51, 0.45, 0.41, 0.36,
                    0.34, 0.30, 0.26, 0.22, 0.17, 0.12, 0.10, 0.10,
                    0.10, 0.10, 0.09, 0.09, 0.10]) 
                    # Holtan værdier for infiltration (growth index) - Sæsonafhængig


Holtan_a = 0.6  # Vegetationsparameter for given afgrød/

# Jordegenskaber

if Tylle == True:
    b = 4.21
    psi_entry = -2.74
    theta_s = 0.435
    theta_wp = 0.069
    theta_airdry = 0.002827


theta_fc = theta_s*(psi_entry/-100)**(1/b)  # Water content at FC [cm^3/cm^3]

if g_dir == -1:
    theta_initial = theta_fc
else:
    theta_initial = theta_hot

if Poulsen == True:
    K_s = 10**(2.8*np.log10(theta_s-theta_fc)+4.3)/(24*60)
theta = np.ones(n_cells)*theta_initial
v_upper = np.max(GI_idle)*Holtan_a *((theta_s - theta[0])*dz)**1.4+K_s
theta_sat = np.ones(n_cells)*theta_s
theta_sat[ror_celle:(ror_celle+10)] = round(1-(rhob/2.65),3)
K_sat = np.ones(n_cells)*K_s
K_sat[ror_celle:(ror_celle+10)] = K_s*((2.65-rhob)/(2.65-rhoref))**3*(rhob/rhoref)**(-2)
K = np.zeros((1, n_cells))[0]
K_L = np.zeros((1, n_cells))[0]
alpha_L = np.zeros((1, n_cells))[0]  # Local slope of the line intersecting [-]
psi = np.zeros((1, n_cells))[0]
v_int = np.zeros((1, n_cells))[0]
alpha_int = np.zeros((1, n_cells))[0]
K_int = np.zeros((1, n_cells))[0]
K_list = np.zeros((t_n//stepsave, n_cells))
v_list = np.zeros((t_n//stepsave, n_cells))
K_L_list = np.zeros((t_n//stepsave, n_cells))
psi_list = np.zeros((t_n//stepsave, n_cells))
theta_list = np.zeros((t_n//stepsave, n_cells))
inflow = np.zeros((t_n))
theta_list[0] = theta
overflow_Q = np.zeros((t_n))
results = np.zeros((t_n//stepsave, n_cells))
d = regn_array[0] # Ponded water at time 0
F = 0 # Accumulated infiltration
I_acc = np.zeros((t_n))
I_acc_theta = np.zeros((t_n))
trtrt = []
ror_celle = 40
A = 2.5*np.pi * 0.08 * 12 / 30
for t in range(1, t_n):
    # Calculation of variables in nodal points
    #index = np.searchsorted(theta_retention, theta, side='left')
    psi = psi_entry*(theta/theta_sat)**(-b)
    K = K_sat * (psi_entry/psi)**(2+3/b)  # K(ψ)
    alpha_L = (2+3/b)/(-psi)
    K_L = K* np.exp(-alpha_L * psi)
    # Calculation of variables at the control surfaces
    for i in range(0, n_cells):
        # Calc values for interfaces between nodal points (control surfaces)
        if i == n_cells-1:
            # Lower boundary condition
            if g_dir == -1:
                v_int[i] = - (K[i]-K_s)/(np.log(K[i]/K_s)) \
                    * (-psi[i]/dz-1) if K[i] != K_s else v_int[i-1]
            elif g_dir == 1:
                v_int[i] = 0 if v_int[i-1] > 0 else v_int[i-1]
                #v_int[i] = - (K[i]-K_s)/(np.log(K[i]/K_s)) \
                 #   * (-psi[i]/dz-1) if K[i] != K_s else v_int[i-1]
            
        else:
            # Slope
            alpha_int[i] = (alpha_L[i] + alpha_L[i+1]) / 2

            # Hydraulic conductivity
            K_int[i] = (K_L[i] + K_L[i+1]) / 2

            # Velocity
            K_n1 = K_int[i] * np.exp(alpha_int[i] * psi[i])
            K_n2 = K_int[i] * np.exp(alpha_int[i] * psi[i+1])
            if g_dir != 0:
                v_int[i] = - (K_n2 - K_n1) / (np.exp(alpha_int[i] * dz) - 1) \
                    + K_int[i] * np.exp(alpha_int[i] * psi[i])
            else:
                # Horizontal flux, Moldrup et al. 1993
                v_int[i] = - (K_n2 - K_n1) / (alpha_int[i] * dz)
    for i in range(0, n_cells):
        # In the first box (i = 0), we need to specify the velocity from the
        # upper boundary condition (v_upper)
        if i == 0:
            v_upper = np.max(GI_idle) * Holtan_a *((theta_s - theta[i])*dz)**1.4+K_s #Holtan infiltration
            tr += regn_array[t]*dt/dz #Nedbør til topreservoir
            if g_dir == -1:
                tr = 0 #Topreservoir er tomt, hvis mod graviation er fokus
            theta[i] = theta[i] - (v_int[i] - min(tr*dz/dt, v_upper)) * dt / dz # Infiltration, mindste led 
            #theta[i] = theta[i] -(v_int[i] - v_upper) * dt/dz
            #theta[i] = theta[i] - (v_int[i] - v_upper)*dt/dz
            F += min(tr*dz/dt, v_upper) * dt / dz
            F_acc.append(F)
            tr -= min(tr*dz/dt, v_upper) * dt / dz
            #tr -= v_upper * dt / dz
            
            tt.append(tr)
            # For the rest of the boxes, we use the calculated velocities at
            # the control surfaces.
        else:
            theta[i] = theta[i] - (v_int[i] - v_int[i-1]) * dt/dz

    while max(theta/theta_sat) > 1: # Tjekker for overmættede vandceller og smider overskydende vand i ovenstående celle
        index_overflow = np.where(theta > theta_sat)[0]
        for overflow in index_overflow:
            overflow_amount = theta[overflow] - theta_sat[overflow]
            theta[overflow] = theta_sat[overflow]
            if overflow < 139:
                v_int[overflow] -= overflow_amount * dz/dt
            if overflow == 0:
                tr += overflow_amount
            else:
                theta[overflow-1] += overflow_amount
    


    kriterie = -5

    if abs(theta[ror_celle]-theta_s) < 10**(kriterie): #Tjekker om cellerne omkring røret er vandmættede
        Pressure_Gradient = (abs((theta[ror_celle::-1]-theta_s)) < 10**(kriterie))
        cum_rvs = Pressure_Gradient.cumsum()
        result = ( cum_rvs - np.maximum.accumulate((~Pressure_Gradient)*cum_rvs ))
        trtrt.append(max(result))
        if (abs(theta[ror_celle:(ror_celle-10):-1] - theta_s) < 10**(kriterie)).all():
            inflow[t] = K_s * dt/dz * A
        elif (abs(theta[ror_celle:(ror_celle-9):-1] - theta_s) < 10**(kriterie)).all():
            inflow[t] = K_s * dt/dz * 0.71 * A
        elif (abs(theta[ror_celle:(ror_celle-8):-1] - theta_s) < 10**(kriterie)).all():
            inflow[t] = K_s * dt/dz * 0.59 * A
        elif (abs(theta[ror_celle:(ror_celle-7):-1] - theta_s) < 10**(kriterie)).all():
            inflow[t] = K_s * dt/dz * 0.49 * A
        elif (abs(theta[ror_celle:(ror_celle-6):-1] - theta_s) < 10**(kriterie)).all():
            inflow[t] = K_s * dt/dz * 0.41 * A
        elif (abs(theta[ror_celle:(ror_celle-5):-1] - theta_s) < 10**(kriterie)).all():
            inflow[t] = K_s * dt/dz * 0.33 * A
        elif (abs(theta[ror_celle:(ror_celle-4):-1] - theta_s) < 10**(kriterie)).all():
            inflow[t] = K_s * dt/dz * 0.26 * A
        elif (abs(theta[ror_celle:(ror_celle-3):-1] - theta_s) < 10**(kriterie)).all():
            inflow[t] = K_s * dt/dz * 0.19 * A
        elif (abs(theta[ror_celle:(ror_celle-2):-1] - theta_s) < 10**(kriterie)).all():
            inflow[t] = K_s * dt/dz * 0.13 * A
        elif abs(theta[ror_celle] - theta_s) < 10**(kriterie):
            inflow[t] = K_s * dt/dz * 0.06 * A
        if max(result) <= ror_celle:
            #theta[10-max(result)] -= K_s * 51.2/500 * dt/dz
            theta[ror_celle-max(result)] -= inflow[t]
        else:
            #tr -= K_s * 51.2/500 * dt/dz
            if inflow[t] > tr:
                Overskud = inflow[t] - tr
                tr = 0
                theta[0] -= Overskud
            else:
                tr -= inflow[t]
        #trtrt.append(result)   
    else:
        Pressure_Gradient = 0
        result = 0
        inflow[t] = 0
      
    theta_list[t//stepsave] = theta


    # Store variables in nested arrays
    K_list[t//stepsave] = K
    v_list[t//stepsave] = v_int
    K_L_list[t//stepsave] = K_L
    psi_list[t//stepsave] = psi

np.savez(f"Model_resultater_{Parameter}{globals()[Parameter]}.npz", Inflow = inflow*dz*500, theta = theta_list, Nedbor = regn_array * dt * 500, Topreservoir = np.array(tt) * dz * 500, Model_vand = theta_list * 500 * dz, v_list = v_list, psi_list = psi_list)