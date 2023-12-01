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
c_c2 = Node(m = m_c_c/2, W = W_c, y0 = T0_c_c2) 
c_m1 = Node(T0 = T0_c_m)
n = Node(y0 = n_frac0)
C1 = Node(y0 = C0[0])
rho = Node(y0 = 0.0)

ARE.addNodes([c_f1,c_f2,c_t1,c_c1,c_c2,c_m1,n,C1,rho])

# define dynamics
c_f1.set_dTdt_BulkFlow(source = (hx_fh1_f2.y(t-tau_hx_c_f)+hx_fh2_f2.y(t-tau_hx_c_f))/2) 
c_f1.set_dTdt_Internal(source = n.y, k = k_f1*P/mcp_f_c)
c_f1.set_dTdt_convective(source = [c_t1.y], hA_mcp = [hA_ft_c/mcp_f_c])

c_f2.set_dTdt_BulkFlow(source = c_f1.y) 
c_f2.set_dTdt_Internal(source = n.y, k = k_f2*P/mcp_f_c)
c_f2.set_dTdt_convective(source = [c_t1.y], hA_mcp = [hA_ft_c/mcp_f_c])

# reactivity 

# dC_i/dt (precursor concentrations)
dC1 = n()*beta[0]/Lam - lam[0]*C1() - C1()/tau_c + C1(t-tau_l)*np.exp(-lam[0]*tau_l)/tau_c                       # C1: y(27)

# reactivity
drho = (a_f/2)*(dc_f1 + dc_f2)+(a_b)*(dc_m)+(a_c/2)*(dc_c1+dc_c2)           # rho: y(33)

# instantiate jitcdde object
DDE = jitcdde([dc_f1,dc_f2,dc_t1,dc_c1,dc_c2,dc_m,dn,dC1,drho])

# set initial conditions
DDE.constant_past([T0_c_f1,T0_c_f2,T0_c_t1,T0_c_c1,T0_c_c2,T0_c_m+50,n_frac0,C0[0],0.0])

T = np.arange(0.0,250,0.01)
sol_jit = []
for t_x in T:
    sol_jit.append(DDE.integrate(t_x))

print(sol_jit[-1][6])
