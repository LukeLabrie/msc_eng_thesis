# %%

from variables_minimal import *
from jitcdde import jitcdde, t
import numpy as np
from heat_transfer import *
from parameters import *
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
        self.dndt = 0.0
        self.dcdt = 0.0
        self.drdt = 0.0
        self.y = None           # temperature (K), to be assigned by System class
    
    def set_dTdt_bulkFlow(self, source):
        '''
        Energy from bulk flow
        source: source node (state variable y(i) or constant)
        dumped: if 'from node' is a constant (indicates dumping instead of 
                recirculation), this needs to be set to true
        '''
        if isinstance(source, float):
                self.BulkFlow = (source-self.y)*self.W/self.m
        elif isinstance(source, type(y)):
                self.BulkFlow = (source()-self.y())*self.W/self.m

    def set_dTdt_internal(self, source: y, k: float):
        '''
        Energy from fission
        n: fractional neutron density 
        P: nominal power (MW)
        mcp: product of node mass (kg) and node specific heat capacity (MJ/(kg*K))
        k: fraction of power generation in node 
        '''
        self.internal = k*source()
    

    def set_dTdt_convective(self, source: list, hA_mcp: list):
        '''
        Energy from convective heat transfer
        a: from node(s) (state variable(s) y(i))
        hA_mcp: ratio of [convective heat transfer coefficient(s) * wetted area(s) (MW/C)]
                to [product of node mass (kg) and node specific heat capcity (MJ)/(kg*k)]
        '''
        for i in range(len(source)):
                self.convective += hA_mcp[i]*(source[i]()-self.y())
    
    def set_dndt(self, r: y, beta_eff: float, Lambda: float, lam: list, C: list):
         '''
         Point kinetics equation for neutron population
         r: reactivity, rho, (state variable y(i))
         beta_eff: delayed neutron fraction
         Lambda: prompt neutron lifetime (s)
         lam: decay constants for associated precursor groups
         C: precursor groups (list of state variables y(i))
         '''
         precursors = 0.0
         for g in enumerate(lam):
              precursors += g[1]*C[g[0]]()
         self.dndt = (r()-beta_eff)*self.y()/Lambda + precursors

    def set_dcdt(self, 
                 n: y, 
                 beta: float, 
                 Lamda: float, 
                 lam: float, 
                 t_c: float, 
                 t_l: float):
         '''
         Precursor concentration
         n: fractional neutron concentration (state variable y(i))
         beta: fraction for precuror group
         Lambda: prompt neutron lifetime 
         lam: decay constant for precursor group
         t_c: core transit time (s)
         t_l: loop transit time (s)
         '''
         source = n()*beta/Lam
         decay = lam*self.y()
         outflow = self.y()/tau_c
         inflow = self.y(t-tau_l)*np.exp(-lam*tau_l)/tau_c
         self.dcdt = source - decay - outflow + inflow 
    
    def set_drdt(self, sources: list, coeffs: list):
         '''
         Reactivity feedback 
         sources: list of derivatives of feedback sources (dy(i)/dt)
         coeffs: list of respective feedback coefficients
         '''
         fb = 0.0
         for s in enumerate(sources):
              fb += s[1]()*coeffs[s[0]]
         self.drdt = fb
    
    def dydt(self):
          y1 = self.BulkFlow
          y2 = self.internal
          y3 = self.convective
          y4 = self.dndt
          y5 = self.dcdt
          y6 = self.drdt
          sum = y1 + y2 + y3 + y4 + y5 + y6
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
                self.nNodes += 1
          return None


# ARE system        
ARE = System()

# instantiate nodes
# core nodes 
c_f1 = Node(m = m_f_c/2, W = W_f, y0 = T0_c_f1)
c_f2 = Node(m = m_f_c/2, W = W_f, y0 = T0_c_f2)
c_t1 = Node(y0 = T0_c_t1)
c_c1 = Node(m = m_c_c/2, W = W_c, y0 = T0_c_c1)
c_c2 = Node(m = m_c_c/2, W = W_c, y0 = T0_c_c2) 
c_m1 = Node(y0 = T0_c_m+50)
n = Node(y0 = n_frac0)
C1 = Node(y0 = C0[0])
rho = Node(y0 = 0.0)

ARE.addNodes([c_f1,c_f2,c_t1,c_c1,c_c2,c_m1,n,C1,rho])

# define dynamics
# core nodes
c_f1.set_dTdt_bulkFlow(source = 500) 
c_f1.set_dTdt_internal(source = n.y, k = k_f1*P/mcp_f_c)
c_f1.set_dTdt_convective(source = [c_t1.y], hA_mcp = [hA_ft_c/mcp_f_c])

c_f2.set_dTdt_bulkFlow(source = c_f1.y) 
c_f2.set_dTdt_internal(source = n.y, k = k_f2*P/mcp_f_c)
c_f2.set_dTdt_convective(source = [c_t1.y], hA_mcp = [hA_ft_c/mcp_f_c])

c_t1.set_dTdt_convective(source = [c_f1.y, c_f2.y, c_c1.y, c_c2.y], 
                         hA_mcp = [hA_ft_c/mcp_t_c, hA_ft_c/mcp_t_c, hA_tc_c/mcp_t_c,hA_tc_c/mcp_t_c])

c_c1.set_dTdt_bulkFlow(source = 500)
c_c1.set_dTdt_convective(source = [c_t1.y], hA_mcp = [hA_tc_c/mcp_c_c])

c_c2.set_dTdt_bulkFlow(source = c_c1.y)
c_c2.set_dTdt_convective(source = [c_t1.y], hA_mcp = [hA_tc_c/mcp_c_c])

c_m1.set_dTdt_internal(source = n.y, k = k_m*P/mcp_m_c)
c_m1.set_dTdt_convective(source = [c_c1.y, c_c2.y], hA_mcp = [hA_mc_c/mcp_m_c]*2)

n.set_dndt(rho.y, beta_t, Lam, [lam[0]], [C1.y])
C1.set_dcdt(n.y, beta[0], Lam, lam[0], tau_c, tau_l)
rho.set_drdt([c_f1.dydt,c_f2.dydt,c_m1.dydt,c_c1.dydt,c_c2.dydt],[a_f/2,a_f/2,a_b,a_c/2,a_c/2])

# instantiate jitcdde object
DDE = jitcdde([n.dydt() for n in ARE.nodes])

# set initial conditions
DDE.constant_past([n.y0 for n in ARE.nodes])



#  %%


T = np.arange(0.0,250,0.01)
sol_jit = []
for t_x in T:
    sol_jit.append(DDE.integrate(t_x))

print(sol_jit[-1][6])
