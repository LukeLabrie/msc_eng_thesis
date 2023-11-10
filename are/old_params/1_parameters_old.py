import numpy as np
import pandas as pd
import math
pi = math.pi

#################################
# ARE
#################################
def iwc_to_Pa(h20):
    '''
    Inches water column -> Pa
    '''
    return h20*248.8

def cfm_to_m3s(v):
    '''
    Cubic ft/min -> m3/s
    '''
    return v/2119

# fuel density (kg/m^3) ORNL-1845 pg. 113
def fuel_density(temp):
    return 1000*3.98-0.00093*temp


# F to K conversion
def F_to_K(tempF):
    return (tempF-32)*5/9+273.15

R = 8.314  # ideal gas constant
def mass_flow(m,P,V,T):
    '''
    m: molar mass of gas
    P: pressure
    V: volume
    T: temperature
    '''
    return m*P*V/(R*T)       # helium mass flow (kg/s)

# thermal feedback (1/Kelvin, temperature provided in Kelvin) ORNL-1845 pg. 115
a_f = (-9.8e-5)*5/9
a_b = (1.1e-5)*5/9

# operating conditions taken from 25-hr Xenon run Exp. H-8
# temperatures 
Vf_c = 1.78/35.315     # fuel volume in core (m^3) ORNL-1845 pg. 108
Vc_c = 10/35.315       # coolant volums in core (m^3) ORNL-1845 pg. 108

T0_c_f1 = F_to_K(1150)     # core fuel inlet temp (K) ORNL-1845 pg. 121
T0_c_f2 = F_to_K(1450)     # core fuel outlet temp (K) could be 1480 if taken from previous page ORNL-1845 pg. 121

T0_c_c1 = F_to_K(1105)     # core coolant inlet temp (K) ORNL-1845 pg. 121
T0_c_c2 = F_to_K(1235)     # core coolant outlet temp (K) ORNL-1845 pg. 121

T0_hfh_f1 = F_to_K(1450)   # fuel-helium hx fuel inlet temp (K) ORNL-1845 pg. 121
T0_hfh_f2 = F_to_K(1150)   # fuel-helium hx fuel oulet temp (K) ORNL-1845 pg. 121

T0_hfh_h1 = F_to_K(180)    # fuel-helium hx helium inlet temp (K) ORNL-1845 pg. 121
T0_hfh_h2 = F_to_K(620)    # fuel-helium hx helium outlet temp (K) ORNL-1845 pg. 121

T0_hhwf_h1 = F_to_K(620)    # helium-water hx (fuel loop) helium inlet temp (K) ORNL-1845 pg. 121
T0_hhwf_h2 = F_to_K(180)    # helium-water hx (fuel loop) helium outlet temp (K) ORNL-1845 pg. 121

T0_hhwf_w1 = F_to_K(70)    # helium-water hx (fuel loop) water inlet temp (K) ORNL-1845 pg. 121 
T0_hhwf_w2 = F_to_K(135)   # helium-water hx (fuel loop) water water outlet temp (K) ORNL-1845 pg. 121

T0_hch_c1 = F_to_K(1235)  # coolant-helium hx coolant inlet temp (K) ORNL-1845 pg. 121 
T0_hch_c2 = F_to_K(1105)  # coolant-helium hx coolant outlet temp (K) ORNL-1845 pg. 121 

T0_hch_h1 = F_to_K(170)  # coolant-helium hx helium inlet temp (K) ORNL-1845 pg. 122 
T0_hch_h2 = F_to_K(1020)  # coolant-helium hx helium outlet temp (K) ORNL-1845 pg. 122 

T0_hhwc_h1 = F_to_K(1020)    # helium-water hx (coolant loop) helium inlet temp (K) ORNL-1845 pg. 122
T0_hhwc_h2 = F_to_K(170)    # helium-water hx (coolant loop) helium outlet temp (K) ORNL-1845 pg. 122

T0_hhwc_w1 = F_to_K(70)    # helium-water hx (coolant loop) water inlet temp (K) ORNL-1845 pg. 121 
T0_hhwc_w2 = F_to_K(100)   # helium-water hx (coolant loop) water water outlet temp (K) ORNL-1845 pg. 121

# flow rates 
F_c_f = 40/15850       # core fuel flow rate (m^3/s) ORNL-1845 pg. 121
F_c_c = 224/15850      # core coolant flow rate (m^3/s) ORNL-1845 pg. 121

F_hfh_h1 = (7300*2)/15850    # fuel-helium hx helium fuel inlet flow rate (m^3/s) ORNL-1845 pg. 121
F_hfh_h2 = (12300*2)/15850    # fuel-helium hx helium fuel outlet flow rate (m^3/s) ORNL-1845 pg. 121

F_hhw_h1 = (12300*2)/15850    # helium-water hx helium fuel inlet flow rate (m^3/s) ORNL-1845 pg. 121
F_hhw_h2 = (7300*2)/15850    # helium-water hx helium fuel outlet flow rate (m^3/s) ORNL-1845 pg. 121

# conversion to mass flow for helium
m_H = 0.004 # molar mass of helium (kg/mol)

# Nominal Power 
P = 2.12              # steady-state power (MW) ORNL-1845 pg. 58

# fuel dimensions
T_fuel_avg = F_to_K(1311)                  # ORNL 1845 pg. 58
V_fuel = 52071.248849/2e6                  # CAD model (m^3)
A_fuel = (69872.856584/2-(6*6.996))/10000 # CAD model (m^2) 

# tube dimensions
V_tubes = 5453.961/1e6                  # CAD model (m^3)
A_tubes = (73445.338-(12*7.728))/10000  # CAD model (m^2) 

# fuel-helium heat exchanger 
L_eff = 93.65/3.281                            # effective tube length (m) ORNL-1535 pg. 47
V_f_hx = (L_eff/1550)*(math.pi*((1.0/2)-0.109)**2)/2 # volume in tubes (m^3) ORNL-1535 pg. 47
m_f_hx = V_f_hx * fuel_density(T_fuel_avg)

# coolant dimensions 
V_coolant = (248933.207)/2e6            # CAD model (m^3)

# fuel mass flow rate 
m_f_c = fuel_density(T_fuel_avg)*V_fuel
W_f = F_c_f * fuel_density(T_fuel_avg)

# coolant mass flow rate
rho_c = 1000*0.78                # coolant density (kg/m^3) 
m_c_c = rho_c*V_coolant      # coolant mass (kg)
W_c = F_c_c * rho_c              # coolant mass flow rate (kg/s)

# delays
tau_hx_c_f = 10.00 # fuel-helium hx to core delay (unknown)
tau_hx_c_c = 10.00 # coolant-helium hx to core delay (unknown)
tau_c_hx_f = 10.00 # coolant-helium hx to core delay (unknown)

# wights
k_f1 = 0.465
k_f2 = 0.465
k_b = 1-(k_f1+k_f2)
k_1 = 0.5
k_2 = 1-k_1

# heat transfer
# fuel -> tubes (core)
h_ft = 6.480e-01  # heat transfer*area coefficient from primary to tubes (MW/C) ORNL-TM-1647 p.3
hA_ft_c = h_ft*A_fuel

# fuel -> tubes (hx)
hA_ft_hx = 5.724*(0.0010550559)*(9/5) # ORNL-1535 pg. 43 CHECK CONVERSION

# tubes -> coolant
h_tc = 3.060E-01 # heat transfer*area coefficient from tubes to secondary (MW/C) ORNL-TM-1647 p.3
hA_tc = h_tc*A_tubes

rho_inconel = 8.5*1000          # inconel density (kg/m^3)
m_t     = V_tubes*rho_inconel   # mass of tubes (kg)  
scp_t = 4.44e-4                 # specific heat capacity of inconel 600 (MJ/kg-C) (https://www.specialmetals.com/documents/technical-bulletins/inconel/inconel-alloy-600.pdf)
mcp_t_c   = m_t*scp_t;          # from ratio of (A_phe/mcp_t)msbr = (A_phe/mcp_t)msre m_tn*cp_tn; % mass*(heat capacity) of tubes per lump in MW-s/°C 

V_t_hx = (L_eff/1550)*math.pi*((1.0/2))**2 - V_f_hx
mcp_t_hx = V_t_hx*rho_inconel

# heat exchanger helium (fuel loop)
rho_h = 0.167                                 # NEEDS TO BE TEMPERATURE DEPENDENT
m_h_hxfh = (((24.75*26.25)/61020)-V_t_hx)*rho_h/2 # hx dimension minus tube volume, factor of two for two nodes
hA_th = (0.955+0.8+0.58)*(0.0010550559)*(9/5)      # ORNL-1535 pg. 43 CHECK CONVERSION

m_h_hxhw = V_f_hx*rho_h

# specific heat fuel
scp_f = 1.9665e-3  # specific heat capacity of fuel salt (MJ/kg-C) ORNL-TM-0728 p.8
mcp_f_c = scp_f*m_f_c
mcp_f_hx = scp_f*m_f_hx

# specific heat coolant
scp_c = 1.256e-3  # specific heat capacity of cooolant (MJ/kg-C) ORNL-1845 p.113
mcp_c_c = scp_c*m_c_c 
m_c_hx = V_f_hx * rho_c      # assume same volume
mcp_h_c = m_c_hx*scp_c

scp_h = 5.1932e-3            # specigic heat capacity of helium (MJ/kg-C)
mcp_h_hxfh = m_h_hxfh*scp_h

# fuel->helium hx
P_hfh_h1 = iwc_to_Pa(2.8)
V_hfh_h1 = cfm_to_m3s(7300)
W_hfh_h1 = mass_flow(m_H,P_hfh_h1,V_hfh_h1,F_to_K(180))

P_hfh_h2 = iwc_to_Pa(1.5)
V_hfh_h2 = cfm_to_m3s(12300)
W_hfh_h2 = mass_flow(m_H,P_hfh_h1,V_hfh_h1,F_to_K(620))

# coolant-> helium hx
P_hch_h1 = iwc_to_Pa(4)
V_hch_h1 = cfm_to_m3s(2000)
W_hch_h1 = mass_flow(m_H,P_hfh_h1,V_hfh_h1,F_to_K(170))

P_hch_h2 = iwc_to_Pa(2)
V_hch_h2 = cfm_to_m3s(4700)
W_hch_h2 = mass_flow(m_H,P_hfh_h1,V_hfh_h1,F_to_K(1020))

# helium->water hx (fuel loop)
P_hhwf_h1 = iwc_to_Pa(1.5)
V_hhwf_h1 = cfm_to_m3s(12300)
W_hhwf_h1 = mass_flow(m_H,P_hfh_h1,V_hfh_h1,F_to_K(620))

P_hhwf_h2 = iwc_to_Pa(1.2)
V_hhwf_h2 = cfm_to_m3s(7300)
W_hhwf_h2 = mass_flow(m_H,P_hfh_h1,V_hfh_h1,F_to_K(180))

W_hhwf_w = 65/15850                # water flow (m^3/s)
m_w = (((24.75*26.25)/61020)-V_t_hx)*998/2 
scp_w = 4181
mcp_w = m_w*scp_w

# helium->water hx (coolant loop)
P_hhwc_h1 = iwc_to_Pa(2.0)
V_hhwc_h1 = cfm_to_m3s(4700)
W_hhwc_h1 = mass_flow(m_H,P_hfh_h1,V_hfh_h1,F_to_K(1020))

P_hhwc_h2 = iwc_to_Pa(0.0)
V_hhwc_h2 = cfm_to_m3s(2000)
W_hhwc_h2 = mass_flow(m_H,P_hfh_h1,V_hfh_h1,F_to_K(170))

W_hhwc_w = 77/15850  

# NEUTRONICS DATA
tau_l = 16.73  # ORNL-TM-0728 %16.44; % (s)
tau_c = 8.46  # ORNL-TM-0728 %8.460; % (s)
n_frac0 = 1.0  # initial fractional neutron density n/n0 (n/cm^3/s)
Lam = 2.400E-04  # mean generation time ORNL-TM-1070 p.15 U235
# Lam = 4.0E-04;  # mean generation time ORNL-TM-1070 p.15 U233
lam = np.array([1.240E-02, 3.05E-02, 1.11E-01, 3.01E-01, 1.140E+00, 3.014E+00])
beta = np.array([0.000223, 0.001457, 0.001307, 0.002628, 0.000766, 0.00023])  # U235
# beta = np.array([0.00023, 0.00079, 0.00067, 0.00073, 0.00013, 0.00009])  # U233
beta_t = np.sum(beta)  # total delayed neutron fraction MSRE
rho_0 = beta_t-sum(np.divide(beta,1+np.divide(1-np.exp(-lam*tau_l),lam*tau_c))) # reactivity change in going from stationary to circulating fuel
C0 = beta / Lam * (1.0 / (lam - (np.exp(-lam * tau_l) - 1.0) / tau_c))

# Tube nodes temps
T0_c_t1 = ((T0_c_f1+T0_c_c1)/2)
T0_hfh_t1 = ((T0_hfh_f1+T0_hfh_h1)/2)
T0_hch_t1 = ((T0_hch_c1+T0_hch_h1)/2)
T0_hhwf_t1 = ((T0_hhwf_h1+T0_hhwf_w1)/2)
T0_hhwc_t1 = ((T0_hhwc_h1+T0_hhwc_w1)/2)

T0_c_b = T_fuel_avg