from heat_transfer import *
from variables_refactor import *
from parameters import *
from jitcdde import y, t
from variables import *
from chspy import CubicHermiteSpline
from jitcdde import y, t, jitcdde_input, input
import matplotlib.pyplot as plt

class Node:
    def __init__(self,
                 m: float = 0.0,
                 W: float = 0.0,
                 T0: float = 0.0) -> None:
        self.m = m              # mass (kg)
        self.W = W              # mass flow rate (kg/s)
        self.T0 = T0            # initial temperature (K)
        self.BulkFlow = 0.0
        self.internal = 0.0
        self.convective = 0.0
        self.T0 = 0.0
        self.T = None           # temperature (K), to be assigned by system class
    
    def set_dTdt_BulkFlow(self, source: y, dumped: bool = False):
        '''
        Energy from bulk flow
        source: source node (state variable y(i))
        dumped: if 'from node' is a constant (indicates dumping instead of 
                recirculation), this needs to be set to true
        '''
        if (dumped):
                self.BulkFlow = (source)*self.W/self.m
        else: 
                self.BulkFlow = (source-self.T)*self.W/self.m
        return None

    def set_dTdt_Internal(self, n: y, P: float, mcp: float, k: float = 1.0):
        '''
        Energy from fission
        n: fractional neutron density 
        P: nominal power (MW)
        mcp: product of node mass (kg) and node specific heat capacity (MJ/(kg*K))
        k: fraction of power generation in node 
        '''
        self.internal = k*P*n/mcp
        return None
    

    def set_dTdt_convective(self, source: list, hA_mcp: list):
        '''
        Energy from convective heat transfer
        a: from node(s) (state variable(s) y(i))
        hA_mcp: ratio of [convective heat transfer coefficient(s) * wetted area(s) (MW/C)]
                to [product of node mass (kg) and node specific heat capcity (MJ)/(kg*k)]
        '''
        for i in range(len(source)):
                self.convective += hA_mcp[i]*(source[i]-self.T)
        return None
    
    def dTdt(self):
          bf = self.BulkFlow
          int = self.internal
          conv = self.convective
          sum = bf + int + conv
          return sum


# core nodes 
c_f1 = Node(m = m_f_c/2, W = W_f, T0 = T0_c_f1)
c_f2 = Node(m = m_f_c/2, W = W_f, T0 = T0_c_f1)


# ARE.addNode(c_f1)

# c_f1.set_dTdt_BulkFlow(source = (hx_fh1_f2(t-tau_hx_c_f)+hx_fh2_f2(t-tau_hx_c_f))/2) 
# c_f1.set_dTdt_Internal(n = n(), P = P, mcp = mcp_f_c, k = k_f1)
# c_f1.set_dTdt_convective(source = [c_t1()], hA_mcp = [hA_ft_c/mcp_f_c])



dc_f2 = dT_bulkFlow(W_f, m_f_c/2, T_c_f1(), c_f2()) + \
        dT_internal(k_f1, P, mcp_f_c, n()) + \
        dT_convective([c_t1()],c_f2(),[hA_ft_c/mcp_f_c])

# tubes
dc_t1= dT_convective([T_c_f1(),c_f2(),c_c1(),c_c2()],c_t1(),[hA_ft_c/mcp_t_c,hA_ft_c/mcp_t_c,hA_tc_c/mcp_t_c,hA_tc_c/mcp_t_c])

# coolant 
dc_c1 = dT_bulkFlow(W_c,m_c_c/2,(hx_ch1_c2(t-tau_hx_c_c)+hx_ch2_c2(t-tau_hx_c_c))/2,c_c1()) + \
        dT_convective([c_t1(),c_m()],c_c1(),[hA_tc_c/mcp_c_c,hA_mc_c/mcp_c_c])

dc_c2 = dT_bulkFlow(W_c,m_c_c/2,c_c1(),c_c2()) + \
        dT_convective([c_t1(),c_m()],c_c2(),[hA_tc_c/mcp_c_c,hA_mc_c/mcp_c_c])

# moderator 
dc_m =  dT_internal(k_m,P,mcp_m_c,n()) + \
        dT_convective([c_c1(),c_c2()],c_m(),[hA_mc_c/mcp_m_c,hA_mc_c/mcp_m_c])

dc_f2 = dT_bulkFlow(W_f, m_f_c/2, T_c_f1(), c_f2()) +                                             dT_internal(k_f1, P, mcp_f_c, n()) + dT_convective([c_t1()],c_f2(),[hA_ft_c/mcp_f_c]) 

# tubes
dc_t1= dT_convective([T_c_f1(),c_f2(),c_c1(),c_c2()],c_t1(),[hA_ft_c/mcp_t_c,hA_ft_c/mcp_t_c,hA_tc_c/mcp_t_c,hA_tc_c/mcp_t_c])

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


# COOLANT-HELIUM HX1
# Fuel Nodes
dhx_ch1_c1 = dT_bulkFlow(W_c, m_c_hx, c_c2(t-tau_c_hx_f), hx_ch1_c1()) + F*dT_convective([hx_ch1_t1()], hx_ch1_c1(), [hA_ct_hx/mcp_c_hxch])
dhx_ch1_c2 = dT_bulkFlow(W_c, m_c_hx, hx_ch1_c1(), hx_ch1_c2()) +        F*dT_convective([hx_ch1_t1()], hx_ch1_c2(), [hA_ct_hx/mcp_c_hxch])

# Tubes
dhx_ch1_t1 = F*dT_convective([hx_ch1_c1(),hx_ch1_c2(),hx_ch1_h1(),hx_ch1_h2()],hx_ch1_t1(),[h/mcp_t_hxch for h in [hA_ct_hx,hA_ct_hx,hA_th_hxch,hA_th_hxch]])

# Helium
dhx_ch1_h1 = dT_bulkFlow(W_h_ch, m_h_hxch/2, hx_hwc1_h2(t-tau_h), hx_ch1_h1()) + F*dT_convective([hx_ch1_t1()], hx_ch1_h1(), [hA_th_hxch/mcp_h_hxch])
dhx_ch1_h2 = dT_bulkFlow(W_h_ch, m_h_hxch/2, hx_ch1_h1(), hx_ch1_h2()) +         F*dT_convective([hx_ch1_t1()], hx_ch1_h2(), [hA_th_hxch/mcp_h_hxch])

# COOLANT-HELIUM HX2
# Fuel Nodes
dhx_ch2_c1 = dT_bulkFlow(W_c, m_c_hx, c_c2(t-tau_c_hx_f), hx_ch2_c1()) + F*dT_convective([hx_ch2_t1()], hx_ch2_c1(), [hA_ct_hx/mcp_c_hxch])
dhx_ch2_c2 = dT_bulkFlow(W_c, m_c_hx, hx_ch2_c1(), hx_ch2_c2()) +        F*dT_convective([hx_ch2_t1()], hx_ch2_c2(), [hA_ct_hx/mcp_c_hxch])

# Tubes
dhx_ch2_t1 = dT_convective([hx_ch2_c1(),hx_ch2_c2(),hx_ch2_h1(),hx_ch2_h2()],hx_ch2_t1(),[h/mcp_t_hxch for h in [hA_ct_hx,hA_ct_hx,hA_th_hxch,hA_th_hxch]])

# Helium
dhx_ch2_h1 = dT_bulkFlow(W_h_ch, m_h_hxch/2, hx_hwc2_h2(t-tau_h), hx_ch2_h1()) + F*dT_convective([hx_ch2_t1()], hx_ch2_h1(), [hA_th_hxch/mcp_h_hxch])
dhx_ch2_h2 = dT_bulkFlow(W_h_ch, m_h_hxch/2, hx_ch2_h1(), hx_ch2_h2()) +         F*dT_convective([hx_ch2_t1()], hx_ch2_h2(), [hA_th_hxch/mcp_h_hxch])


# HELIUM-WATER HX1 (FUEL LOOP)
# Helium
dhx_hwf1_h1 = dT_bulkFlow(W_h_fh, m_h_hxhw/2, hx_fh1_h2(), hx_hwf1_h1()) +  F*dT_convective([hx_hwf1_t1()], hx_hwf1_h1(), [hA_ht_hxhw/mcp_h_hxhw])
dhx_hwf1_h2 = dT_bulkFlow(W_h_fh, m_h_hxhw/2, hx_hwf1_h1(), hx_hwf1_h2()) + F*dT_convective([hx_hwf1_t1()], hx_hwf1_h2(), [hA_ht_hxhw/mcp_h_hxhw])

# Tubes
dhx_hwf1_t1 = F*dT_convective([hx_hwf1_h1(),hx_hwf1_h2(),hx_hwf1_w1(),hx_hwf1_w2()],hx_hwf1_t1(),[h/mcp_t_hxhw for h in [hA_ht_hxhw,hA_ht_hxhw,hA_tw_hxhw,hA_tw_hxhw]])

# Water
dhx_hwf1_w1 = dT_bulkFlow(W_hhwf_w, m_w/2, T0_hhwf_w1-hx_hwf1_w1(), hx_hwf1_w1(),dumped=True) + F*dT_convective([hx_hwf1_t1()], hx_hwf1_w1(), [hA_tw_hxhw/mcp_w])
dhx_hwf1_w2 = dT_bulkFlow(W_hhwf_w, m_w/2, hx_hwf1_w1(), hx_hwf1_w2()) +                        F*dT_convective([hx_hwf1_t1()], hx_hwf1_w2(), [hA_tw_hxhw/mcp_w])

# HELIUM-WATER HX2 (FUEL LOOP)
# Helium
dhx_hwf2_h1 = dT_bulkFlow(W_h_fh, m_h_hxhw/2, hx_fh2_h2(), hx_hwf2_h1()) +  F*dT_convective([hx_hwf2_t1()], hx_hwf2_h1(), [hA_ht_hxhw/mcp_h_hxhw])
dhx_hwf2_h2 = dT_bulkFlow(W_h_fh, m_h_hxhw/2, hx_hwf2_h1(), hx_hwf2_h2()) + F*dT_convective([hx_hwf2_t1()], hx_hwf2_h2(), [hA_ht_hxhw/mcp_h_hxhw])

# Tubes
dhx_hwf2_t1 = F*dT_convective([hx_hwf2_h1(),hx_hwf2_h2(),hx_hwf2_w1(),hx_hwf2_w2()],hx_hwf2_t1(),[h/mcp_t_hxhw for h in [hA_ht_hxhw,hA_ht_hxhw,hA_tw_hxhw,hA_tw_hxhw]])

# Water
dhx_hwf2_w1 = dT_bulkFlow(W_hhwf_w, m_w/2, T0_hhwf_w1-hx_hwf2_w1(), hx_hwf2_w1(),dumped=True) + F*dT_convective([hx_hwf2_t1()], hx_hwf2_w1(), [hA_tw_hxhw/mcp_w])
dhx_hwf2_w2 = dT_bulkFlow(W_hhwf_w, m_w/2, hx_hwf2_w1(), hx_hwf2_w2()) +                        F*dT_convective([hx_hwf2_t1()], hx_hwf2_w2(), [hA_tw_hxhw/mcp_w])


# HELIUM-WATER HX1 (COOLANT LOOP)
# Helium
dhx_hwc1_h1 = dT_bulkFlow(W_h_ch, m_h_hxhwc/2, hx_ch1_h2(), hx_hwc1_h1()) + F*dT_convective([hx_hwc1_t1()], hx_hwc1_h1(), [hA_ht_hxhwc/mcp_h_hxhwc])
dhx_hwc1_h2 = dT_bulkFlow(W_h_ch, m_h_hxhwc/2, hx_hwc1_h1(), hx_hwc1_h2()) + F*dT_convective([hx_hwc1_t1()], hx_hwc1_h2(), [hA_ht_hxhwc/mcp_h_hxhwc])

# Tubes
dhx_hwc1_t1 = F*dT_convective([hx_hwc1_h1(),hx_hwc1_h2(),hx_hwc1_w1(),hx_hwc1_w2()],hx_hwc1_t1(),[h/mcp_t_hxhwc for h in [hA_ht_hxhwc,hA_ht_hxhwc,hA_tw_hxhwc,hA_tw_hxhwc]])

# Water
dhx_hwc1_w1 = dT_bulkFlow(W_hhwc_w, m_w_hxhwc/2, T0_hhwf_w1-hx_hwc1_w1(), hx_hwc1_w1(), dumped=True) + F*dT_convective([hx_hwc1_t1()], hx_hwc1_w1(), [hA_tw_hxhwc/mcp_w_hxhwc])
dhx_hwc1_w2 = dT_bulkFlow(W_hhwc_w, m_w_hxhwc/2, hx_hwc1_w1(), hx_hwc1_w2()) + F*dT_convective([hx_hwc1_t1()], hx_hwc1_w2(), [hA_tw_hxhwc/mcp_w_hxhwc])


# HELIUM-WATER HX2 (COOLANT LOOP)
# Helium
dhx_hwc2_h1 = dT_bulkFlow(W_h_ch, m_h_hxhwc/2, hx_ch2_h2(), hx_hwc2_h1()) + F*dT_convective([hx_hwc2_t1()], hx_hwc2_h1(), [hA_ht_hxhwc/mcp_h_hxhwc])
dhx_hwc2_h2 = dT_bulkFlow(W_h_ch, m_h_hxhwc/2, hx_hwc2_h1(), hx_hwc2_h2()) + F*dT_convective([hx_hwc2_t1()], hx_hwc2_h2(), [hA_ht_hxhwc/mcp_h_hxhwc])

# Tubes
dhx_hwc2_t1 = F*dT_convective([hx_hwc2_h1(),hx_hwc2_h2(),hx_hwc2_w1(),hx_hwc2_w2()],hx_hwc2_t1(),[h/mcp_t_hxhwc for h in [hA_ht_hxhwc,hA_ht_hxhwc,hA_tw_hxhwc,hA_tw_hxhwc]])

# Water
dhx_hwc2_w1 = dT_bulkFlow(W_hhwc_w, m_w_hxhwc/2, T0_hhwf_w1-hx_hwc2_w1(), hx_hwc2_w1(), dumped=True) + F*dT_convective([hx_hwc2_t1()], hx_hwc2_w1(), [hA_tw_hxhwc/mcp_w_hxhwc])
dhx_hwc2_w2 = dT_bulkFlow(W_hhwc_w, m_w_hxhwc/2, hx_hwc2_w1(), hx_hwc2_w2()) + F*dT_convective([hx_hwc2_t1()], hx_hwc2_w2(), [hA_tw_hxhwc/mcp_w_hxhwc])

# external input
def rho_insert(t):
    if (t<t_ins):
        return 0.0
    elif (t<(t_ins+insert_duration)):
        return ((t-t_ins))*(inserted/insert_duration) # linear
    else:
        return inserted

rho_spline = CubicHermiteSpline(n=1)
rho_spline.from_function(rho_insert, times_of_interest = T)
rho_ext = input(0)


dn = ((rho()+rho_ext)-beta_t)*n()/Lam+lam[0]*C1()+lam[1]*C2()+lam[2]*C3()+lam[3]*C4()+lam[4]*C5()+lam[5]*C6()           # n (no source insertion): n()

# dC_i/dt (precursor concentrations)
dC1 = n()*beta[0]/Lam - lam[0]*C1() - C1()/tau_c + C1(t-tau_l)*np.exp(-lam[0]*tau_l)/tau_c                       # C1: y(27)
dC2 = n()*beta[1]/Lam - lam[1]*C2() - C2()/tau_c + C2(t-tau_l)*np.exp(-lam[1]*tau_l)/tau_c                       # C2: y(28)
dC3 = n()*beta[2]/Lam - lam[2]*C3() - C3()/tau_c + C3(t-tau_l)*np.exp(-lam[2]*tau_l)/tau_c                       # C3: y(29)
dC4 = n()*beta[3]/Lam - lam[3]*C4() - C4()/tau_c + C4(t-tau_l)*np.exp(-lam[3]*tau_l)/tau_c                       # C4: y(30)
dC5 = n()*beta[4]/Lam - lam[4]*C5() - C5()/tau_c + C5(t-tau_l)*np.exp(-lam[4]*tau_l)/tau_c                       # C5: y(31)
dC6 = n()*beta[5]/Lam - lam[5]*C6() - C6()/tau_c + C6(t-tau_l)*np.exp(-lam[5]*tau_l)/tau_c                       # C6: y(32)

# reactivity 
drho = (a_f/2)*(c_f1.dTdt() + dc_f2)+(a_b)*(dc_m)+(a_c/2)*(dc_c1+dc_c2)           # rho: y(33)

# Initial values & solve


# instantiate jitcdde object
DDE = jitcdde_input([c_f1.dTdt(),dc_f2,dc_t1,dc_c1,dc_c2,dc_m,
                     dhx_fh1_f1,dhx_fh1_f2,dhx_fh1_t1,dhx_fh1_h1,dhx_fh1_h2,
                     dhx_fh2_f1,dhx_fh2_f2,dhx_fh2_t1,dhx_fh2_h1,dhx_fh2_h2,
                     dhx_ch1_c1,dhx_ch1_c2,dhx_ch1_t1,dhx_ch1_h1,dhx_ch1_h2,
                     dhx_ch2_c1,dhx_ch2_c2,dhx_ch2_t1,dhx_ch2_h1,dhx_ch2_h2,
                     dhx_hwf1_h1,dhx_hwf1_h2,dhx_hwf1_t1,dhx_hwf1_w1,dhx_hwf1_w2,
                     dhx_hwf2_h1,dhx_hwf2_h2,dhx_hwf2_t1,dhx_hwf2_w1,dhx_hwf2_w2,
                     dhx_hwc1_h1,dhx_hwc1_h2,dhx_hwc1_t1,dhx_hwc1_w1,dhx_hwc1_w2,
                     dhx_hwc2_h1,dhx_hwc2_h2,dhx_hwc2_t1,dhx_hwc2_w1,dhx_hwc2_w2,
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


#DDE.set_integration_parameters(atol=1e-50, rtol=1e-20, first_step=1.0, min_step=1e-11, max_step=10.0, decrease_threshold=1.1, 
#                           increase_threshold=0.5, safety_factor=0.9, max_factor=5.0, min_factor=0.2, pws_factor=3, 
#                           pws_atol=0.0, pws_rtol=1e-05, pws_max_iterations=10, pws_base_increase_chance=0.1, pws_fuzzy_increase=False)
#
# DDE.step_on_discontinuities()
sol_jit = []
for t_x in T:
    sol_jit.append(DDE.integrate(t_x))


plt.plot(T,[s[46] for s in sol_jit])
plt.savefig('test_run.png')
print(sol_jit[-1][46])