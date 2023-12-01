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
                 y0: float = 0.0) -> None:
        self.m = m              # mass (kg)
        self.W = W              # mass flow rate (kg/s)
        self.y0 = y0            # initial temperature (K)
        self.BulkFlow = 0.0
        self.internal = 0.0
        self.convective = 0.0
        self.y0 = 0.0
        self.y = None           # temperature (K), to be assigned by System class
    
    def set_dTdt_BulkFlow(self, source: y, dumped: bool = False):
        '''
        Energy from bulk flow
        source: source node (state variable y(i))
        dumped: if 'from node' is a constant (indicates dumping instead of 
                recirculation), this needs to be set to true
        '''
        if (dumped):
                self.BulkFlow = (source.y)*self.W/self.m
        else: 
                self.BulkFlow = (source.y-self.y)*self.W/self.m
        return None

    def set_dTdt_Internal(self, source: y, k: float):
        '''
        Energy from fission
        n: fractional neutron density 
        P: nominal power (MW)
        mcp: product of node mass (kg) and node specific heat capacity (MJ/(kg*K))
        k: fraction of power generation in node 
        '''
        self.internal = k*source.y
        return None
    

    def set_dTdt_convective(self, source: list, hA_mcp: list):
        '''
        Energy from convective heat transfer
        a: from node(s) (state variable(s) y(i))
        hA_mcp: ratio of [convective heat transfer coefficient(s) * wetted area(s) (MW/C)]
                to [product of node mass (kg) and node specific heat capcity (MJ)/(kg*k)]
        '''
        for i in range(len(source)):
                self.convective += hA_mcp[i]*(source[i].y-self.y)
        return None
    
    def dTdt(self):
          T1 = self.BulkFlow
          T2 = self.internal
          T3 = self.convective
          sum = T1 + T2 + T3
          return sum


class System:
     def __init__(self) -> None:
          self.nNodes = 0
          self.nodes = []

     def addNodes(self, newNodes: list):
          index = self.nNodes
          for n in newNodes: 
                n.y = lambda tau=None, i = index: y(i, tau) if tau is not None else y(i)
                self.nodes.append(n)
                nNodes += 1
          return None


# ARE system        
ARE = System()

# instantiate nodes
# core nodes 
c_f1 = Node(m = m_f_c/2, W = W_f, y0 = T0_c_f1)
c_f2 = Node(m = m_f_c/2, W = W_f, y0 = T0_c_f1)
c_t1 = Node(T0 = T0_c_t1)
c_c1 = Node(m = m_c_c/2, W = W_c, y0 = T0_c_c1)
c_c1 = Node(m = m_c_c/2, W = W_c, y0 = T0_c_c2) 
c_m1 = Node(T0 = T0_c_m)

# neutron kinetics 
n = Node(y0 = n_frac0)
C1 = Node(y0 = C0[0])
C2 = Node(y0 = C0[1])
C3 = Node(y0 = C0[2])
C4 = Node(y0 = C0[3])
C5 = Node(y0 = C0[4])
C6 = Node(y0 = C0[5])
rho = Node(y0 = 0.0)

# hx fuel->helium 1
hx_fh1_f1 = Node(m = m_f_hx, W = W_f, y0 = T0_hfh_f1)
hx_fh1_f1 = Node(m = m_f_hx, W = W_f, y0 = T0_hfh_f2)
hx_fh1_t1 = Node(m = m_t_hxfh, y0 = T0_hfh_t1)
hx_fh1_h1 = Node(m = m_h_hxfh/2, W = W_h_fh, y0 = T0_hfh_h1)
hx_fh1_h2 = Node(m = m_h_hxfh/2, W = W_h_fh, y0 = T0_hfh_h2) 

# hx fuel->helium 2
hx_fh2_f1 = Node(m = m_f_hx, W = W_f, y0 = T0_hfh_f1)
hx_fh2_f1 = Node(m = m_f_hx, W = W_f, y0 = T0_hfh_f2)
hx_fh2_t1 = Node(m = m_t_hxfh, y0 = T0_hfh_t1)
hx_fh2_h1 = Node(m = m_h_hxfh/2, W = W_h_fh, y0 = T0_hfh_h1)
hx_fh2_h2 = Node(m = m_h_hxfh/2, W = W_h_fh, y0 = T0_hfh_h2) 

# hx coolant->helium 1
hx_ch1_c1 = Node(m = m_c_hx, W = W_c, y0 = T0_hch_c1)
hx_ch1_c2 = Node(m = m_c_hx, W = W_c, y0 = T0_hch_c2)
hx_ch1_t1 = Node(m = m_t_hxch, y0 = T0_hch_t1)
hx_ch1_h1 = Node(m = m_h_hxch/2, W = W_c, y0 = T0_hch_h1)
hx_ch1_h2 = Node(m = m_h_hxch/2, W = W_c, y0 = T0_hch_h2)

# hx coolant->helium 2
hx_ch2_c1 = Node(m = m_c_hx, W = W_c, y0 = T0_hch_c1)
hx_ch2_c2 = Node(m = m_c_hx, W = W_c, y0 = T0_hch_c2)
hx_ch2_t1 = Node(m = m_t_hxch, y0 = T0_hch_t1)
hx_ch2_h1 = Node(m = m_h_hxch/2, W = W_c, y0 = T0_hch_h1)
hx_ch2_h2 = Node(m = m_h_hxch/2, W = W_c, y0 = T0_hch_h2)

# hx helium->water 1, fuel loop
hx_hwf1_h1 = Node(m = m_h_hxhw/2, W = W_h_fh, y0 = T0_hhwf_h1)
hx_hwf1_h2 = Node(m = m_h_hxhw/2, W = W_h_fh, y0 = T0_hhwf_h2)
hx_hwf1_t1 = Node(m = m_t_hxhw, y0 = T0_hhwf_t1)
hx_hwf1_w1 = Node(m = m_w_hxhw/2, W = W_hhwf_w, y0 = T0_hhwf_w1)
hx_hwf1_w2 = Node(m = m_w_hxhw/2, W = W_hhwf_w, y0 = T0_hhwf_w2)

# hx helium->water 2, fuel loop
hx_hwf1_h1 = Node(m = m_h_hxhw/2, W = W_h_fh, y0 = T0_hhwf_h1)
hx_hwf1_h2 = Node(m = m_h_hxhw/2, W = W_h_fh, y0 = T0_hhwf_h2)
hx_hwf1_t1 = Node(m = m_t_hxhw, y0 = T0_hhwf_t1)
hx_hwf1_w1 = Node(m = m_w_hxhw/2, W = W_hhwf_w, y0 = T0_hhwf_w1)
hx_hwf1_w2 = Node(m = m_w_hxhw/2, W = W_hhwf_w, y0 = T0_hhwf_w2)

# hx helium->water 1, coolant loop
hx_hwc1_h1 = Node(m = m_h_hxhwc, W = W_h_ch, y0 = T0_hhwc_h1)
hx_hwc1_h2 = Node(m = m_h_hxhwc, W = W_h_ch, y0 = T0_hhwc_h1)
hx_hwc1_t1 = Node(m = m_h_hxhwc, y0 = T0_hhwc_t1)
hx_hwc1_w1 = Node(m = m_w_hchwc, W = W_hhwc_w, y0 = T0_hhwc_w1)
hx_hwc1_w2 = Node(m = m_w_hchwc, W = W_hhwc_w, y0 = T0_hhwc_w2)

# hx helium->water 2, coolant loop
hx_hwc1_h1 = Node(m = m_h_hxhwc, W = W_h_ch, y0 = T0_hhwc_h1)
hx_hwc1_h2 = Node(m = m_h_hxhwc, W = W_h_ch, y0 = T0_hhwc_h1)
hx_hwc1_t1 = Node(m = m_h_hxhwc, y0 = T0_hhwc_t1)
hx_hwc1_w1 = Node(m = m_w_hchwc, W = W_hhwc_w, y0 = T0_hhwc_w1)
hx_hwc1_w2 = Node(m = m_w_hchwc, W = W_hhwc_w, y0 = T0_hhwc_w2)

# add nodes to system
ARE.addNodes([c_f1,dc_f2,dc_t1,dc_c1,dc_c2,dc_m,
              dn,dC1,dC2,dC3,dC4,dC5,dC6,drho,
              dhx_fh1_f1,dhx_fh1_f2,dhx_fh1_t1,dhx_fh1_h1,dhx_fh1_h2,
              dhx_fh2_f1,dhx_fh2_f2,dhx_fh2_t1,dhx_fh2_h1,dhx_fh2_h2,
              dhx_ch1_c1,dhx_ch1_c2,dhx_ch1_t1,dhx_ch1_h1,dhx_ch1_h2,
              dhx_ch2_c1,dhx_ch2_c2,dhx_ch2_t1,dhx_ch2_h1,dhx_ch2_h2,
              dhx_hwf1_h1,dhx_hwf1_h2,dhx_hwf1_t1,dhx_hwf1_w1,dhx_hwf1_w2,
              dhx_hwf2_h1,dhx_hwf2_h2,dhx_hwf2_t1,dhx_hwf2_w1,dhx_hwf2_w2,
              dhx_hwc1_h1,dhx_hwc1_h2,dhx_hwc1_t1,dhx_hwc1_w1,dhx_hwc1_w2,
              dhx_hwc2_h1,dhx_hwc2_h2,dhx_hwc2_t1,dhx_hwc2_w1,dhx_hwc2_w2,
                ])

# define dynamics
c_f1.set_dTdt_BulkFlow(source = (hx_fh1_f2.y(t-tau_hx_c_f)+hx_fh2_f2.y(t-tau_hx_c_f))/2) 
c_f1.set_dTdt_Internal(source = n.y, k = k_f1*P/mcp_f_c)
c_f1.set_dTdt_convective(source = [c_t1.y], hA_mcp = [hA_ft_c/mcp_f_c])

c_f2.set_dTdt_BulkFlow(source = c_f1.y) 
c_f2.set_dTdt_Internal(source = n.y, k = k_f2*P/mcp_f_c)
c_f2.set_dTdt_convective(source = [c_t1.y], hA_mcp = [hA_ft_c/mcp_f_c])

# try with a minimal example