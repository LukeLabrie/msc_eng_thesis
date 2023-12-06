# %%

from variables_minimal import *
from jitcdde import jitcdde, t
import numpy as np
from heat_transfer import *
from parameters import *
import matplotlib.pyplot as plt

# CORE
# fuel nodes
dc_f1 = dT_bulkFlow(W_f, m_f_c/2, 500, c_f1()) + dT_internal(k_f1, P, mcp_f_c, n()) + dT_convective([c_t1()],c_f1(),[hA_ft_c/mcp_f_c]) 
dc_f2 = dT_bulkFlow(W_f, m_f_c/2, c_f1(), c_f2()) +                                             dT_internal(k_f1, P, mcp_f_c, n()) + dT_convective([c_t1()],c_f2(),[hA_ft_c/mcp_f_c]) 

# tubes
dc_t1= dT_convective([c_f1(),c_f2(),c_c1(),c_c2()],c_t1(),[hA_ft_c/mcp_t_c,hA_ft_c/mcp_t_c,hA_tc_c/mcp_t_c,hA_tc_c/mcp_t_c])

# coolant 
dc_c1 = dT_bulkFlow(W_c,m_c_c/2, 500, c_c1()) + dT_convective([c_t1(),c_m()],c_c1(),[hA_tc_c/mcp_c_c,hA_mc_c/mcp_c_c])
dc_c2 = dT_bulkFlow(W_c,m_c_c/2,c_c1(),c_c2()) + dT_convective([c_t1(),c_m()],c_c2(),[hA_tc_c/mcp_c_c,hA_mc_c/mcp_c_c]) 

# moderator 
dc_m =  dT_internal(k_m,P,mcp_m_c,n()) + F*dT_convective([c_c1(),c_c2()],c_m(),[hA_mc_c/mcp_m_c,hA_mc_c/mcp_m_c])

dn = ((rho())-beta_t)*n()/Lam+lam[0]*C1()

# reactivity 

# dC_i/dt (precursor concentrations)
dC1 = n()*beta[0]/Lam - lam[0]*C1() - C1()/tau_c + C1(t-tau_l)*np.exp(-lam[0]*tau_l)/tau_c                       # C1: y(27)

# reactivity
drho = (a_f/2)*(dc_f1 + dc_f2)+(a_b)*(dc_m)+(a_c/2)*(dc_c1+dc_c2)           # rho: y(33)

# instantiate jitcdde object
DDE = jitcdde([dc_f1,dc_f2,dc_t1,dc_c1,dc_c2,dc_m,dn,dC1,drho])

# set initial conditions
DDE.constant_past([T0_c_f1,T0_c_f2,T0_c_t1,T0_c_c1,T0_c_c2,T0_c_m+50,n_frac0,C0[0],0.0])

DDE.get_state()
# %%

T = np.arange(0.0,250,0.01)
sol_jit = []
for t_x in T:
    sol_jit.append(DDE.integrate(t_x))

plt.plot(T,[s[6] for s in sol_jit])
plt.show()
print(sol_jit[-1][6])
