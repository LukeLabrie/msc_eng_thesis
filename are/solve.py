from heat_transfer import *
from variables import *
from parameters import *
from jitcdde import y, t

class Node:
    def __init__(self,
                 T: y,
                 m: float = 0.0,
                 W: float = 0.0) -> None:
        self.T = T              # temperature (K)
        self.m = m              # mass (kg)
        self.W = W              # mass flow rate (kg/s)
        self.BulkFlow = 0.0
        self.internal = 0.0
        self.convective = 0.0
    
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

# system for solver
dydt = []
y0 = []

# instantiate nodes
# core 
c_f1 = Node(T = c_f1(), m = m_f_c/2, W = W_f)
c_f1.set_dTdt_BulkFlow(source = (hx_fh1_f2(t-tau_hx_c_f)+hx_fh2_f2(t-tau_hx_c_f))/2) 
c_f1.set_dTdt_Internal(n = n(), P = P, mcp = mcp_f_c, k = k_f1)
c_f1.set_dTdt_convective(source = [c_t1()], hA_mcp = [hA_ft_c/mcp_f_c])
dydt.append(c_f1.dTdt())
y0.append(T0_c_f1)

 #.... continue   


# instantiate jitcdde object
#DDE = jitcdde_input([dc_f1,dc_f2,dc_t1,dc_c1,dc_c2,dc_m,
#                     dhx_fh1_f1,dhx_fh1_f2,dhx_fh1_t1,dhx_fh1_h1,dhx_fh1_h2,
#                     dhx_fh2_f1,dhx_fh2_f2,dhx_fh2_t1,dhx_fh2_h1,dhx_fh2_h2,
#                     dhx_ch1_c1,dhx_ch1_c2,dhx_ch1_t1,dhx_ch1_h1,dhx_ch1_h2,
#                     dhx_ch2_c1,dhx_ch2_c2,dhx_ch2_t1,dhx_ch2_h1,dhx_ch2_h2,
#                     dhx_hwf1_h1,dhx_hwf1_h2,dhx_hwf1_t1,dhx_hwf1_w1,dhx_hwf1_w2,
#                     dhx_hwf2_h1,dhx_hwf2_h2,dhx_hwf2_t1,dhx_hwf2_w1,dhx_hwf2_w2,
#                     dhx_hwc1_h1,dhx_hwc1_h2,dhx_hwc1_t1,dhx_hwc1_w1,dhx_hwc1_w2,
#                     dhx_hwc2_h1,dhx_hwc2_h2,dhx_hwc2_t1,dhx_hwc2_w1,dhx_hwc2_w2,
#                     dn,dC1,dC2,dC3,dC4,dC5,dC6,drho],
#                     rho_spline)

# set initial conditions
#DDE.constant_past([T0_c_f1,T0_c_f2,T0_c_t1,T0_c_c1,T0_c_c2,T0_c_m+50,
 #                  T0_hfh_f1,T0_hfh_f2,T0_hfh_t1,T0_hfh_h1,T0_hfh_h2,
 #                  T0_hfh_f1,T0_hfh_f2,T0_hfh_t1,T0_hfh_h1,T0_hfh_h2,
 #                  T0_hch_c1,T0_hch_c2,T0_hch_t1,T0_hch_h1,T0_hch_h2,
 #                  T0_hch_c1,T0_hch_c2,T0_hch_t1,T0_hch_h1,T0_hch_h2,
 #                  T0_hhwf_h1,T0_hhwf_h2,T0_hhwf_t1,T0_hhwf_w1,T0_hhwf_w2,
#                   T0_hhwf_h1,T0_hhwf_h2,T0_hhwf_t1,T0_hhwf_w1,T0_hhwf_w2,
#                   T0_hhwc_h1,T0_hhwc_h2,T0_hhwc_t1,T0_hhwc_w1,T0_hhwc_w2,
#                   T0_hhwc_h1,T0_hhwc_h2,T0_hhwc_t1,T0_hhwc_w1,T0_hhwc_w2,
#                   n_frac0,C0[0],C0[1],C0[2],C0[3],C0[4],C0[5],0.0
#                    ])



# define system


# CORE
# fuel nodes
dc_f1 = dT_bulkFlow(W_f, m_f_c/2, (hx_fh1_f2(t-tau_hx_c_f)+hx_fh2_f2(t-tau_hx_c_f))/2,c_f1()) +  \
        dT_internal(k_f1, P, mcp_f_c, n()) + \
        dT_convective([c_t1()],c_f1(),[hA_ft_c/mcp_f_c])

dc_f2 = dT_bulkFlow(W_f, m_f_c/2, c_f1(), c_f2()) + \
        dT_internal(k_f1, P, mcp_f_c, n()) + \
        dT_convective([c_t1()],c_f2(),[hA_ft_c/mcp_f_c])

# tubes
dc_t1= dT_convective([c_f1(),c_f2(),c_c1(),c_c2()],c_t1(),[hA_ft_c/mcp_t_c,hA_ft_c/mcp_t_c,hA_tc_c/mcp_t_c,hA_tc_c/mcp_t_c])

# coolant 
dc_c1 = dT_bulkFlow(W_c,m_c_c/2,(hx_ch1_c2(t-tau_hx_c_c)+hx_ch2_c2(t-tau_hx_c_c))/2,c_c1()) + \
        dT_convective([c_t1(),c_m()],c_c1(),[hA_tc_c/mcp_c_c,hA_mc_c/mcp_c_c])

dc_c2 = dT_bulkFlow(W_c,m_c_c/2,c_c1(),c_c2()) + \
        dT_convective([c_t1(),c_m()],c_c2(),[hA_tc_c/mcp_c_c,hA_mc_c/mcp_c_c])

# moderator 
dc_m =  dT_internal(k_m,P,mcp_m_c,n()) + \
        dT_convective([c_c1(),c_c2()],c_m(),[hA_mc_c/mcp_m_c,hA_mc_c/mcp_m_c])

