from jitcdde import jitcdde, y, t, jitcdde_input, input

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
