# %%
from jitcdde import jitcdde, y, t, jitcdde_input, input
from parameters_expanded import *
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
from variables import *
from chspy import CubicHermiteSpline

# %%
# heat transfer functions
def dT_convective(a: list, b: y, hA_mcp: list):
    '''
    Energy from convective heat transfer
    a: from node(s) (state variable(s) y(i))
    b: to node (state variable y(i))
    hA_mcp: ratio of [convective heat transfer coefficient(s) * wetted area(s) (MW/C)]
            to [product of node mass (kg) and node specific heat capcity (MJ)/(kg*k)]
    f: fractional adjustments
    '''
    tot = 0.0
    for i in range(len(a)):
        tot += hA_mcp[i]*(a[i]-b)
    return tot

def dT_bulkFlow(W: float, m: float, a: y, b: y, dumped: bool = False):
    '''
    Energy from bulk flow
    W: mass flow rate (kg/s)
    m: node mass (kg)
    a: Q from node (state variable y(i))
    b: Q to node (state variable y(i))
    dumped: if 'from node' is a constant (indicates dumping instead of 
            recirculation), this needs to be set to true
    '''
    if (dumped):
        return (a)*W/m
    else: 
        return (a-b)*W/m

def dT_internal(k: float, P: float, mcp: float, n: y):
    '''
    Energy from fission
    k: fraction of power generation in node 
    P: nominal power (MW)
    mcp: product of node mass (kg) and node specific heat capacity (MJ/(kg*K))
    n: fractional neutron density 
    '''
    return k*P*n/mcp


# %% [markdown]
# Define system

# %%
# CORE
# fuel nodes
dc_f1 = dT_bulkFlow(W_f, m_f_c/2, (hx_fh1_f2(t-tau_hx_c_f)+hx_fh2_f2(t-tau_hx_c_f))/2,c_f1()) + dT_internal(k_f1, P, mcp_f_c, n()) + F*dT_convective([c_t1()],c_f1(),[hA_ft_c/mcp_f_c]) - arbitrary_removal*c_f1()
dc_f2 = dT_bulkFlow(W_f, m_f_c/2, c_f1(), c_f2()) +                                             dT_internal(k_f1, P, mcp_f_c, n()) + dT_convective([c_t1()],c_f2(),[hA_ft_c/mcp_f_c]) - arbitrary_removal*c_f2()

# tubes
dc_t1= F*dT_convective([c_f1(),c_f2(),c_c1(),c_c2()],c_t1(),[hA_ft_c/mcp_t_c,hA_ft_c/mcp_t_c,hA_tc_c/mcp_t_c,hA_tc_c/mcp_t_c]) - arbitrary_removal*c_t1()

# coolant 
dc_c1 = dT_bulkFlow(W_c,m_c_c/2,(hx_ch1_c2(t-tau_hx_c_c)+hx_ch2_c2(t-tau_hx_c_c))/2,c_c1()) + F*dT_convective([c_t1(),c_m()],c_c1(),[hA_tc_c/mcp_c_c,hA_mc_c/mcp_c_c]) - arbitrary_removal*c_c1()
dc_c2 = dT_bulkFlow(W_c,m_c_c/2,c_c1(),c_c2()) +                                              F*dT_convective([c_t1(),c_m()],c_c2(),[hA_tc_c/mcp_c_c,hA_mc_c/mcp_c_c]) - arbitrary_removal*c_c2()

# moderator 
dc_m =  dT_internal(k_m,P,mcp_m_c,n()) + F*dT_convective([c_c1(),c_c2()],c_m(),[hA_mc_c/mcp_m_c,hA_mc_c/mcp_m_c])


# %%
# FUEL-HELIUM HX1
dhx_fh1_f1 = dT_bulkFlow(W_f,m_f_hx,c_f2(t-tau_c_hx_f),hx_fh1_f1()) + F*dT_convective([hx_fh1_t1()],hx_fh1_f1(),[hA_ft_hx/mcp_f_hx]) 
dhx_fh1_f2 = dT_bulkFlow(W_f,m_f_hx,hx_fh1_f1(),hx_fh1_f2()) +        F*dT_convective([hx_fh1_t1()],hx_fh1_f2(),[hA_ft_hx/mcp_f_hx]) 

# Tubes 
dhx_fh1_t1 = dT_convective([hx_fh1_f1(),hx_fh1_f2(),hx_fh1_h1(),hx_fh1_h2()],hx_fh1_t1(),[h/mcp_t_hx for h in [hA_ft_hx,hA_ft_hx,hA_ht_hx,hA_ht_hx]])

# Helium
dhx_fh1_h1 = dT_bulkFlow(W_h_fh,m_h_hxfh,hx_hwf2_h2(t-tau_h),hx_fh1_h1()) + F*dT_convective([hx_fh1_t1()],hx_fh1_h1(),[hA_ht_hx/mcp_h_hxfh])
dhx_fh1_h2 = dT_bulkFlow(W_h_fh,m_h_hxfh,hx_fh1_h1(),hx_fh1_h2()) +         F*dT_convective([hx_fh1_t1()],hx_fh1_h2(),[hA_ht_hx/mcp_h_hxfh])

# FUEL-HELIUM HX2 Fuel Nodes
dhx_fh2_f1 = dT_bulkFlow(W_f,m_f_hx,c_f2(t-tau_c_hx_f), hx_fh2_f1()) + F*dT_convective([hx_fh2_t1()],hx_fh2_f1(),[hA_ft_hx/mcp_f_hx])
dhx_fh2_f2 = dT_bulkFlow(W_f,m_f_hx,hx_fh2_f1(),hx_fh2_f2()) +         F*dT_convective([hx_fh2_t1()],hx_fh2_f2(),[hA_ft_hx/mcp_f_hx])

# Tubes for FUEL-HELIUM HX2
dhx_fh2_t1 = dT_convective([hx_fh2_f1(),hx_fh2_f2(),hx_fh2_h1(),hx_fh2_h2()],hx_fh2_t1(),[h/mcp_t_hx for h in [hA_ft_hx,hA_ft_hx,hA_ht_hx,hA_ht_hx]])

# Helium for FUEL-HELIUM HX2
dhx_fh2_h1 = dT_bulkFlow(W_h_fh,m_h_hxfh,hx_hwf1_h2(),hx_fh2_h1()) + F*dT_convective([hx_fh2_t1()],hx_fh2_h1(),[hA_ht_hx/mcp_h_hxfh])
dhx_fh2_h2 = dT_bulkFlow(W_h_fh,m_h_hxfh,hx_fh2_h1(),hx_fh2_h2()) +  F*dT_convective([hx_fh2_t1()],hx_fh2_h2(),[hA_ht_hx/mcp_h_hxfh])



# %%
# COOLANT-HELIUM HX1
# Fuel Nodes
T_hch_c1 = dT_bulkFlow(W_c, m_c_hx, c_c2(t-tau_c_hx_f), hx_ch1_c1()) + F*dT_convective([hx_ch1_t1()], hx_ch1_c1(), [hA_ct_hx/mcp_c_hxch])
T_hch_c2 = dT_bulkFlow(W_c, m_c_hx, hx_ch1_c1(), hx_ch1_c2()) +        F*dT_convective([hx_ch1_t1()], hx_ch1_c2(), [hA_ct_hx/mcp_c_hxch])

# Tubes
T_hch_t1 = F*dT_convective([hx_ch1_c1(),hx_ch1_c2(),hx_ch1_h1(),hx_ch1_h2()],hx_ch1_t1(),[h/mcp_t_hxch for h in [hA_ct_hx,hA_ct_hx,hA_th_hxch,hA_th_hxch]])

# Helium
T_hch_h1 = dT_bulkFlow(W_h_ch, m_h_hxch/2, hx_hwc1_h2(t-tau_h), hx_ch1_h1()) + F*dT_convective([hx_ch1_t1()], hx_ch1_h1(), [hA_th_hxch/mcp_h_hxch])
T_hch_h2 = dT_bulkFlow(W_h_ch, m_h_hxch/2, hx_ch1_h1(), hx_ch1_h2()) +         F*dT_convective([hx_ch1_t1()], hx_ch1_h2(), [hA_th_hxch/mcp_h_hxch])

# COOLANT-HELIUM HX2
# Fuel Nodes
T_hch2_c1 = dT_bulkFlow(W_c, m_c_hx, c_c2(t-tau_c_hx_f), hx_ch2_c1()) + F*dT_convective([hx_ch2_t1()], hx_ch2_c1(), [hA_ct_hx/mcp_c_hxch])
T_hch2_c2 = dT_bulkFlow(W_c, m_c_hx, hx_ch2_c1(), hx_ch2_c2()) +        F*dT_convective([hx_ch2_t1()], hx_ch2_c2(), [hA_ct_hx/mcp_c_hxch])

# Tubes
T_hch2_t1 = dT_convective([hx_ch2_c1(),hx_ch2_c2(),hx_ch2_h1(),hx_ch2_h2()],hx_ch2_t1(),[h/mcp_t_hxch for h in [hA_ct_hx,hA_ct_hx,hA_th_hxch,hA_th_hxch]])

# Helium
T_hch2_h1 = dT_bulkFlow(W_h_ch, m_h_hxch/2, hx_hwc2_h2(t-tau_h), hx_ch2_h1()) + F*dT_convective([hx_ch2_t1()], hx_ch2_h1(), [hA_th_hxch/mcp_h_hxch])
T_hch2_h2 = dT_bulkFlow(W_h_ch, m_h_hxch/2, hx_ch2_h1(), hx_ch2_h2()) +         F*dT_convective([hx_ch2_t1()], hx_ch2_h2(), [hA_th_hxch/mcp_h_hxch])


# %%
# HELIUM-WATER HX1 (FUEL LOOP)
# Helium
T_hhwf_h1 = dT_bulkFlow(W_h_fh, m_h_hxhw/2, hx_fh1_h2(), hx_hwf1_h1()) +  F*dT_convective([hx_hwf1_t1()], hx_hwf1_h1(), [hA_ht_hxhw/mcp_h_hxhw])
T_hhwf_h2 = dT_bulkFlow(W_h_fh, m_h_hxhw/2, hx_hwf1_h1(), hx_hwf1_h2()) + F*dT_convective([hx_hwf1_t1()], hx_hwf1_h2(), [hA_ht_hxhw/mcp_h_hxhw])

# Tubes
T_hhwf_t1 = F*dT_convective([hx_hwf1_h1(),hx_hwf1_h2(),hx_hwf1_w1(),hx_hwf1_w2()],hx_hwf1_t1(),[h/mcp_t_hxhw for h in [hA_ht_hxhw,hA_ht_hxhw,hA_tw_hxhw,hA_tw_hxhw]])

# Water
T_hhwf_w1 = dT_bulkFlow(W_hhwf_w, m_w/2, T0_hhwf_w1-hx_hwf1_w1(), hx_hwf1_w1(),dumped=True) + F*dT_convective([hx_hwf1_t1()], hx_hwf1_w1(), [hA_tw_hxhw/mcp_w])
T_hhwf_w2 = dT_bulkFlow(W_hhwf_w, m_w/2, hx_hwf1_w1(), hx_hwf1_w2()) +                        F*dT_convective([hx_hwf1_t1()], hx_hwf1_w2(), [hA_tw_hxhw/mcp_w])

# HELIUM-WATER HX2 (FUEL LOOP)
# Helium
T_hhwf2_h1 = dT_bulkFlow(W_h_fh, m_h_hxhw/2, hx_fh2_h2(), hx_hwf2_h1()) +  F*dT_convective([hx_hwf2_t1()], hx_hwf2_h1(), [hA_ht_hxhw/mcp_h_hxhw])
T_hhwf2_h2 = dT_bulkFlow(W_h_fh, m_h_hxhw/2, hx_hwf2_h1(), hx_hwf2_h2()) + F*dT_convective([hx_hwf2_t1()], hx_hwf2_h2(), [hA_ht_hxhw/mcp_h_hxhw])

# Tubes
T_hhwf2_t1 = F*dT_convective([hx_hwf2_h1(),hx_hwf2_h2(),hx_hwf2_w1(),hx_hwf2_w2()],hx_hwf2_t1(),[h/mcp_t_hxhw for h in [hA_ht_hxhw,hA_ht_hxhw,hA_tw_hxhw,hA_tw_hxhw]])

# Water
T_hhwf2_w1 = dT_bulkFlow(W_hhwf_w, m_w/2, T0_hhwf_w1-hx_hwf2_w1(), hx_hwf2_w1(),dumped=True) + F*dT_convective([hx_hwf2_t1()], hx_hwf2_w1(), [hA_tw_hxhw/mcp_w])
T_hhwf2_w2 = dT_bulkFlow(W_hhwf_w, m_w/2, hx_hwf2_w1(), hx_hwf2_w2()) +                        F*dT_convective([hx_hwf2_t1()], hx_hwf2_w2(), [hA_tw_hxhw/mcp_w])


# %%
# HELIUM-WATER HX1 (COOLANT LOOP)
# Helium
T_hhwc_h1 = dT_bulkFlow(W_h_ch, m_h_hxhwc/2, hx_ch1_h2(), hx_hwc1_h1()) + F*dT_convective([hx_hwc1_t1()], hx_hwc1_h1(), [hA_ht_hxhwc/mcp_h_hxhwc])
T_hhwc_h2 = dT_bulkFlow(W_h_ch, m_h_hxhwc/2, hx_hwc1_h1(), hx_hwc1_h2()) + F*dT_convective([hx_hwc1_t1()], hx_hwc1_h2(), [hA_ht_hxhwc/mcp_h_hxhwc])

# Tubes
T_hhwc_t1 = F*dT_convective([hx_hwc1_h1(),hx_hwc1_h2(),hx_hwc1_w1(),hx_hwc1_w2()],hx_hwc1_t1(),[h/mcp_t_hxhwc for h in [hA_ht_hxhwc,hA_ht_hxhwc,hA_tw_hxhwc,hA_tw_hxhwc]])

# Water
T_hhwc_w1 = dT_bulkFlow(W_hhwc_w, m_w_hxhwc/2, T0_hhwf_w1-hx_hwc1_w1(), hx_hwc1_w1(), dumped=True) + F*dT_convective([hx_hwc1_t1()], hx_hwc1_w1(), [hA_tw_hxhwc/mcp_w_hxhwc])
T_hhwc_w2 = dT_bulkFlow(W_hhwc_w, m_w_hxhwc/2, hx_hwc1_w1(), hx_hwc1_w2()) + F*dT_convective([hx_hwc1_t1()], hx_hwc1_w2(), [hA_tw_hxhwc/mcp_w_hxhwc])


# HELIUM-WATER HX2 (COOLANT LOOP)
# Helium
T_hhwc2_h1 = dT_bulkFlow(W_h_ch, m_h_hxhwc/2, hx_ch2_h2(), hx_hwc2_h1()) + F*dT_convective([hx_hwc2_t1()], hx_hwc2_h1(), [hA_ht_hxhwc/mcp_h_hxhwc])
T_hhwc2_h2 = dT_bulkFlow(W_h_ch, m_h_hxhwc/2, hx_hwc2_h1(), hx_hwc2_h2()) + F*dT_convective([hx_hwc2_t1()], hx_hwc2_h2(), [hA_ht_hxhwc/mcp_h_hxhwc])

# Tubes
T_hhwc2_t1 = F*dT_convective([hx_hwc2_h1(),hx_hwc2_h2(),hx_hwc2_w1(),hx_hwc2_w2()],hx_hwc2_t1(),[h/mcp_t_hxhwc for h in [hA_ht_hxhwc,hA_ht_hxhwc,hA_tw_hxhwc,hA_tw_hxhwc]])

# Water
T_hhwc2_w1 = dT_bulkFlow(W_hhwc_w, m_w_hxhwc/2, T0_hhwf_w1-hx_hwc2_w1(), hx_hwc2_w1(), dumped=True) + F*dT_convective([hx_hwc2_t1()], hx_hwc2_w1(), [hA_tw_hxhwc/mcp_w_hxhwc])
T_hhwc2_w2 = dT_bulkFlow(W_hhwc_w, m_w_hxhwc/2, hx_hwc2_w1(), hx_hwc2_w2()) + F*dT_convective([hx_hwc2_t1()], hx_hwc2_w2(), [hA_tw_hxhwc/mcp_w_hxhwc])



# %%
# reactivity insertion at t = 1000
rho_0 = 0.0027582318523071635+beta_t
t0 = 500.00
tf = 1000.00
T = np.arange(0.0,tf,0.01)
inserted = 2*(4e-3)
def rho_insert(t):
    insert_duration = 0.4/0.011
    if (t<t0):
        return 0.0
    elif (t<(t0+insert_duration)):
        return ((t-t0))*(inserted/insert_duration)
    else:
        return inserted

rho_spline = CubicHermiteSpline(n=1)
rho_spline.from_function(rho_insert, times_of_interest = T)

rho_ext = input(0)

# %%

dn = ((rho()+rho_ext)-beta_t)*n()/Lam+lam[0]*C1()+lam[1]*C2()+lam[2]*C3()+lam[3]*C4()+lam[4]*C5()+lam[5]*C6()           # n (no source insertion): n()

# dC_i/dt (precursor concentrations)
dC1 = n()*beta[0]/Lam - lam[0]*C1() - C1()/tau_c + C1(t-tau_l)*np.exp(-lam[0]*tau_l)/tau_c                       # C1: y(27)
dC2 = n()*beta[1]/Lam - lam[1]*C2() - C2()/tau_c + C2(t-tau_l)*np.exp(-lam[1]*tau_l)/tau_c                       # C2: y(28)
dC3 = n()*beta[2]/Lam - lam[2]*C3() - C3()/tau_c + C3(t-tau_l)*np.exp(-lam[2]*tau_l)/tau_c                       # C3: y(29)
dC4 = n()*beta[3]/Lam - lam[3]*C4() - C4()/tau_c + C4(t-tau_l)*np.exp(-lam[3]*tau_l)/tau_c                       # C4: y(30)
dC5 = n()*beta[4]/Lam - lam[4]*C5() - C5()/tau_c + C5(t-tau_l)*np.exp(-lam[4]*tau_l)/tau_c                       # C5: y(31)
dC6 = n()*beta[5]/Lam - lam[5]*C6() - C6()/tau_c + C6(t-tau_l)*np.exp(-lam[5]*tau_l)/tau_c                       # C6: y(32)

# reactivity 
drho = (a_f/2)*(dc_f1 + dc_f2)+(a_b)*(dc_m)+(a_c/2)*(dc_c1+dc_c2)           # rho: y(33)

# %% [markdown]
# Initial values & solve

# %%
# instantiate jitcdde object
DDE = jitcdde_input([dc_f1,dc_f2,dc_t1,dc_c1,dc_c2,dc_m,
                     dhx_fh1_f1,dhx_fh1_f2,dhx_fh1_t1,dhx_fh1_h1,dhx_fh1_h2,
                     dhx_fh2_f1,dhx_fh2_f2,dhx_fh2_t1,dhx_fh2_h1,dhx_fh2_h2,
                     T_hch_c1,T_hch_c2,T_hch_t1,T_hch_h1,T_hch_h2,
                     T_hch2_c1,T_hch2_c2,T_hch2_t1,T_hch2_h1,T_hch2_h2,
                     T_hhwf_h1,T_hhwf_h2,T_hhwf_t1,T_hhwf_w1,T_hhwf_w2,
                     T_hhwf2_h1,T_hhwf2_h2,T_hhwf2_t1,T_hhwf2_w1,T_hhwf2_w2,
                     T_hhwc_h1,T_hhwc_h2,T_hhwc_t1,T_hhwc_w1,T_hhwc_w2,
                     T_hhwc2_h1,T_hhwc2_h2,T_hhwc2_t1,T_hhwc2_w1,T_hhwc2_w2,
                     dn,dC1,dC2,dC3,dC4,dC5,dC6,drho],
                     rho_spline)

# set initial conditions
DDE.constant_past([T0_c_f1,T0_c_f2,T0_c_t1,T0_c_c1,T0_c_c2,T0_c_m+50,
                   T0_hfh_f1,T0_hfh_f2,T0_hfh_t1,T0_hfh_h1,T0_hfh_h2,
                   T0_hfh_f1,T0_hfh_f2,T0_hfh_t1,T0_hfh_h1,T0_hfh_h2,
                   T0_hch_c1,T0_hch_c2,T0_hch_t1,T0_hch_h1,T0_hch_h2,
                   T0_hch_c1,T0_hch_c2,T0_hch_t1,T0_hch_h1,T0_hch_h2,
                   T0_hhwf_h1,T0_hhwf_h2,T0_hhwf_t1,T0_hhwf_w1,T0_hhwf_w2,
                   T0_hhwf_h1,T0_hhwf_h2,T0_hhwf_t1,T0_hhwf_w1,T0_hhwf_w2,
                   T0_hhwc_h1,T0_hhwc_h2,T0_hhwc_t1,T0_hhwc_w1,T0_hhwc_w2,
                   T0_hhwc_h1,T0_hhwc_h2,T0_hhwc_t1,T0_hhwc_w1,T0_hhwc_w2,
                   n_frac0,C0[0],C0[1],C0[2],C0[3],C0[4],C0[5],0.0
                    ])

# %%


#DDE.set_integration_parameters(atol=1e-10, rtol=1e-05, first_step=1.0, min_step=1e-11, max_step=10.0, decrease_threshold=1.1, 
#                           increase_threshold=0.5, safety_factor=0.9, max_factor=5.0, min_factor=0.2, pws_factor=3, 
#                           pws_atol=0.0, pws_rtol=1e-05, pws_max_iterations=10, pws_base_increase_chance=0.1, pws_fuzzy_increase=False)
#
# DDE.step_on_discontinuities()

sol_jit = []
for t_x in T:
    sol_jit.append(DDE.integrate(t_x))

# %%
df_reversed = pd.read_csv("./data/dollar_insert.csv",header=None)
df = df_reversed.iloc[::-1]
df = df.reset_index(drop=True)

# adjust to reported power
#d = df[1][0]-P
#df[1] = [p-d for p in df[1]]

#plt.plot(df[0],df[1])
#plt.xlabel("Time")
#plt.ylabel("Power (MW)")
#plt.title("ARE: 1$ Reactivity Insertion")

# %%
t_before_data = (1110-df[0][0])*60
duration_data = (df.iloc[-1][0]-df[0][0])*60

t_end_data = df.iloc[-1][0]
t_before_sim = t0-t_before_data

T_insert = [t for t in T if (t > (t_before_sim)) and (t < (t_before_sim)+(duration_data))]
i_insert = [t[0] for t in enumerate(T) if (t[1] > (t_before_sim)) and (t[1] < (t_before_sim)+(duration_data))]

adj = (df[0][0])*60-T_insert[0]
T_insert = [t+adj for t in T_insert]

plt.plot(df[0]*60,df[1],label='ORNL-1845')
plt.plot(T_insert,[s[46]*P for s in sol_jit][i_insert[0]:(i_insert[-1]+1)],label=r"simulated") #$\alpha_f'$ = $\frac{\alpha_f}{3.5}$")
plt.legend()
plt.ylabel("Power (MW)")
plt.xlabel("Time(s)")
plt.title("ORNL-1845 vs Simulation")
