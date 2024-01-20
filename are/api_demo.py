from parameters import *
from msrDynamics.objects import Node, System

T_f_in =  Node(m = m_f_in,   # mass
               c = scp_f,    # specific heat capacity 
               y0 = T0_f_in, # initial value
               W = W_f, )    # mass flow
T_f_out = Node(m = m_f_out, c = scp_f, y0 = T0_f_out, W = W_f, )
T_c_in =  Node(m = m_c_in,  c = scp_c, y0 = T0_c_in,  W = W_c)
T_c_out = Node(m = m_c_out, c = scp_c, y0 = T0_c_out, W = W_c,) 
T_t =     Node(m = m_t,     c = scp_t, y0 = T0_t)

# instantiate System instance and add nodes
core = System()
core.add_nodes([T_f_in,T_f_out,T_c_in,T_c_out,T_t])


# fuel bulk flow
T_f_in.set_dTdt_bulkFlow(source = T_f_hx.y(t-tau_hx_c)) 
T_f_out.set_dTdt_bulkFlow(source = T_f_in.y()) 

# fuel convective heat transfer 
T_f_in.set_dTdt_convective(source = [T_t.y()], hA = hA_ft)
T_f_out.set_dTdt_convective(source = [T_t.y()], hA = hA_ft_c)

# tubes convective heat transfer
T_t.set_dTdt_convective(source = [T_f_in.y(), T_f_out.y(), T_c_in.y(), T_c_out.y()], 
                        hA = [hA_ft,hA_ft,hA_tc,hA_tc])

# coolant bulk flow
T_c_in.set_dTdt_bulkFlow(source = T_c_hx.y(t-tau_hx_c))
T_c_out.set_dTdt_bulkFlow(source = T_c_in.y())

# coolant convective heat transfer
T_c_in.set_dTdt_convective(source = T_t.y(), hA = hA_tc)
T_c_out.set_dTdt_convective(source = T_t.y(), hA = hA_tc)