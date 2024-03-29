# %%
from relaxations import *
from parameters import *
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
from scipy.interpolate import CubicSpline

# %%
# read experimental data
df_reversed = pd.read_csv("./data/insertion.csv",header=None)
df = df_reversed.iloc[::-1]
df = df.reset_index(drop=True)

# get indicies for comparison
t_before_data = (1110-df[0][0])*60
duration_data = (df.iloc[-1][0]-df[0][0])*60
t_end_data = df.iloc[-1][0]
t_before_sim = t_ins-t_before_data
T_insert = [t for t in T if (t > (t_before_sim)) and (t < (t_before_sim)+(duration_data))]
i_insert = [t[0] for t in enumerate(T) if (t[1] > (t_before_sim)) and (t[1] < (t_before_sim)+(duration_data))]

adj = (df[0][0])*60-T_insert[0]
df[0] = [(t*60)-adj for t in df[0]]

# adjust to reported power
d = df[1][0]-P
df[1] = [p-d for p in df[1]]

# Set up interpolation
# Assuming df[0] is time and df[1] is the data you want to interpolate
spline = CubicSpline(df[0], df[1])  # Multiplying df[0] by 60 if it's in minutes

# Use the spline to interpolate at the desired times
interpolated_values = spline(T[i_insert[0]:i_insert[-1]+1])


def sumSq_delays(params):
    try:
        # run
        sol_jit = relax_delays(params) # Implement this function to run your model
        
        # calculate error
        simulation_output = [s[46]*P for s in sol_jit][i_insert[0]:(i_insert[-1]+1)]
        error = sum((simulation_output - interpolated_values)**2)  # Sum of squared errors
        return error
    except:
        return float('inf')
    
def sumSq_ht(params):
    try:
        sol_jit = relax_hA_oop(params) # Implement this function to run your model
        
        # calculate error
        T_eq = 200
        i_eq = [i for i in range(len(T)) if T[i] >= T_eq]
        simulation_output = np.array([s[6]*P for s in sol_jit])
        simulation_output = simulation_output[i_eq]
        error = sum((simulation_output - P)**2)  # Sum of squared errors
        return error
    except:
        return float('inf')
    
def sumSq_ht_all(params):
    try:
        sol_jit = relax_hA_oop(params) 

        T = np.arange(0,10000,0.1)
        # calculate error
        T_eq = 9000
        i_eq = [i for i in range(len(T)) if T[i] >= T_eq]

        pow = np.array([s[6]*P for s in sol_jit])[i_eq]
        pow_error = sum((pow - P)**2)  

        fuel_in = np.array([s[15] for s in sol_jit])[i_eq]
        fuel_in_error = sum((fuel_in-F_to_K(1209))**2)

        fuel_out = np.array([s[1] for s in sol_jit])[i_eq]
        fuel_out_error = sum((fuel_out-F_to_K(1522))**2)

        coolant_in = np.array([s[25] for s in sol_jit])[i_eq]
        coolant_in_error = sum((coolant_in-F_to_K(1226))**2)

        coolant_out = np.array([s[4] for s in sol_jit])[i_eq]
        coolant_out_error = sum((coolant_out-F_to_K(1335))**2)

        tot = pow_error + fuel_in_error + fuel_out_error + coolant_in_error + coolant_out_error

        return tot
    except:
        return float('inf')

def sumSq_all(params):
    try:
        sol_jit = relax_all(params) # Implement this function to run your model
        
        # calculate error
        simulation_output = [s[6]*P for s in sol_jit][i_insert[0]:(i_insert[-1]+1)]
        error = sum((simulation_output - interpolated_values)**2)  # Sum of squared errors
        return error
    except:
        print("failed")
        return float('inf')
    

def sumSq_initial(params):
    try:
        # run
        sol_jit = relax_fib_oop(params)
        
        # calculate error
        simulation_output = [s[6]*P for s in sol_jit][i_insert[0]:(i_insert[-1]+1)]
        error = sum((simulation_output - interpolated_values)**2)  # Sum of squared errors
        return error
    except:
        return float('inf')
    
# %%

from scipy.optimize import minimize

def estimate_delays():

    # set bounds
    initial_hx_c_f = 20.0
    hx_c_f_bounds = (10.0, 60.0)
    
    initial_hx_c_c = 20.0
    hx_c_c_bounds = (10.0, 60.0)

    initial_c_hx_f = 20.0
    c_hx_f_bounds = (10.0, 60.0)

    initial_h_loop_f = 1.0
    h_loop_f_bounds = (0.0, 5.0)

    initial_h_btw_f = 0.0
    h_btw_f_bounds = (0.0, 10.0)

    initial_c_hx_c = 20.0
    c_hx_c_bounds = (10.0, 60.0)

    initial_h_loop_c = 1.0
    h_loop_c_bounds = (0.0, 5.0)

    initial_h_btw_c = 0.0
    h_btw_c_bounds = (0.0, 3.0)

    initial_guess = [initial_hx_c_f,initial_hx_c_c,initial_c_hx_f,initial_h_loop_f,
                     initial_h_btw_f,initial_c_hx_c,initial_h_loop_c,initial_h_btw_c]


    bounds = [hx_c_f_bounds,hx_c_c_bounds,c_hx_f_bounds,h_loop_f_bounds,h_btw_f_bounds,
              c_hx_c_bounds,h_loop_c_bounds,h_btw_c_bounds]

    # minimize
    result = minimize(sumSq_delays, initial_guess, bounds=bounds)

    return result

def estimate_ht():

    range_factor = 10

    ft_c = hA_ft_c
    ft_c_lims = (ft_c/range_factor,range_factor*ft_c)

    tc_c = hA_tc_c
    tc_c_lims = (tc_c/range_factor,range_factor*tc_c)

    mc_c = hA_mc_c
    mc_c_lims = (mc_c/range_factor,range_factor*mc_c)

    ft_hx = hA_ft_hx
    ft_hx_lims = (ft_hx/range_factor,range_factor*ft_hx)

    ht_hx = hA_ht_hx
    ht_hx_lims = (ht_hx/range_factor,range_factor*ht_hx)

    ct_hx = hA_ct_hx
    ct_hx_lims = (ct_hx/range_factor,range_factor*ct_hx)

    th_hxch = hA_th_hxch
    th_hxch_lims = (th_hxch/range_factor,range_factor*th_hxch)

    ht_hxhw = hA_ht_hxhw
    ht_hxhw_lims = (ht_hxhw/range_factor,range_factor*ht_hxhw)

    tw_hxhw = hA_tw_hxhw
    tw_hxhw_lims = (tw_hxhw/range_factor,range_factor*tw_hxhw)

    ht_hxhwc = hA_ht_hxhwc
    ht_hxhwc_lims = (ht_hxhwc/range_factor,range_factor*ht_hxhwc)

    tw_hxhwc = hA_tw_hxhwc
    tw_hxhwc_lims = (tw_hxhwc/range_factor,range_factor*tw_hxhwc)

    T0_m = T0_c_m
    T0_m_lims = (T0_c_m-600,T0_c_m+600)


    initial_guess = [ft_c,tc_c,mc_c,ft_hx,ht_hx,ct_hx,th_hxch,ht_hxhw,
                     tw_hxhw,ht_hxhwc,tw_hxhwc,T0_m]
    
    bounds = [ft_c_lims,tc_c_lims,mc_c_lims,ft_hx_lims,ht_hx_lims,ct_hx_lims,
              th_hxch_lims,ht_hxhw_lims,tw_hxhw_lims,ht_hxhwc_lims,tw_hxhwc_lims,
              T0_m_lims]


    result = minimize(
        sumSq_ht_all, 
        initial_guess, 
        bounds=bounds,
        method='Nelder-Mead'
    )

    return result


def estimate_all():

    # set bounds
    a_f0 = -9.8e-5
    a_f_bounds = (-20e-5, 20e-5)

    a_b0 = -1.1e-5
    a_b_bounds = (-5e-5, 5e-5)

    a_c0 = -5.88e-5
    a_c_bounds = (-20e-5, 20e-5)

    ins0 = 4e-3
    ins_bounds = (1e-3, 20e-3)

    b0 = beta_t
    b_bounds = (1e-3,10e-3)

    initial_hx_c_f = 20.0
    hx_c_f_bounds = (10.0, 60.0)
    
    initial_hx_c_c = 20.0
    hx_c_c_bounds = (10.0, 60.0)

    initial_c_hx_f = 20.0
    c_hx_f_bounds = (10.0, 60.0)

    initial_h_loop_f = 1.0
    h_loop_f_bounds = (0.0, 5.0)

    initial_h_btw_f = 0.0
    h_btw_f_bounds = (0.0, 10.0)

    initial_c_hx_c = 20.0
    c_hx_c_bounds = (10.0, 60.0)

    initial_h_loop_c = 1.0
    h_loop_c_bounds = (0.0, 5.0)

    initial_h_btw_c = 0.0
    h_btw_c_bounds = (0.0, 3.0)

    ft_c = 0.025249076460017884
    ft_c_lims = (ft_c/10,10*ft_c)

    tc_c = 0.0146062610110109899
    tc_c_lims = (tc_c/10,10*tc_c)

    mc_c = 0.00092151047071557199
    mc_c_lims = (mc_c/10,10*mc_c)

    ft_hx = 0.006092568792077965
    ft_hx_lims = (ft_hx/10,10*ft_hx)

    ht_hx = 0.0014320505785117184
    ht_hx_lims = (ht_hx/10,10*ht_hx)

    ct_hx = 0.0101010207925026710110
    ct_hx_lims = (ct_hx/10,10*ct_hx)

    th_hxch = 0.0004489850066827337
    th_hxch_lims = (th_hxch/10,10*th_hxch)

    ht_hxhw = 0.004725554058974901
    ht_hxhw_lims = (ht_hxhw/10,10*ht_hxhw)

    tw_hxhw = 0.3439054124906395
    tw_hxhw_lims = (tw_hxhw/10,10*tw_hxhw)

    ht_hxhwc = 0.0004752963985070788
    ht_hxhwc_lims = (ht_hxhwc/10,10*ht_hxhwc)

    tw_hxhwc = 0.0893816147929607
    tw_hxhwc_lims = (tw_hxhwc/10,10*tw_hxhwc)

    Lam_guess = Lam
    Lam_lims = (Lam/3,3*Lam)

    l1_guess = lam[0]
    l1_lims = (l1_guess/5,5*l1_guess)

    l2_guess = lam[1]
    l2_lims = (l2_guess/5,5*l2_guess)

    l3_guess = lam[2]
    l3_lims = (l3_guess/5,5*l3_guess)

    l4_guess = lam[3]
    l4_lims = (l4_guess/5,5*l4_guess)

    l5_guess = lam[4]
    l5_lims = (l5_guess/5,5*l5_guess)

    l6_guess = lam[5]
    l6_lims = (l6_guess/5,5*l6_guess)

    initial_guess = [a_f0,a_b0,a_c0,ins0,b0,initial_hx_c_f,initial_hx_c_c,
                     initial_c_hx_f,initial_h_loop_f,initial_h_btw_f,
                     initial_c_hx_c,initial_h_loop_c,initial_h_btw_c,ft_c,
                     tc_c,mc_c,ft_hx,ht_hx,ct_hx,th_hxch,ht_hxhw,tw_hxhw,
                     ht_hxhwc,tw_hxhwc,Lam_guess,l1_guess,l2_guess,l3_guess,
                     l4_guess,l5_guess,l6_guess]
    
    bounds = [a_f_bounds,a_b_bounds,a_c_bounds,ins_bounds,b_bounds,
              hx_c_f_bounds,hx_c_c_bounds,c_hx_f_bounds,h_loop_f_bounds,
              h_btw_f_bounds,c_hx_c_bounds,h_loop_c_bounds,h_btw_c_bounds, 
              ft_c_lims,tc_c_lims, mc_c_lims,ft_hx_lims,ht_hx_lims,ct_hx_lims,
              th_hxch_lims,ht_hxhw_lims,tw_hxhw_lims,ht_hxhwc_lims,tw_hxhwc_lims,
              Lam_lims,l1_lims,l2_lims,l3_lims,l4_lims,l5_lims,l6_lims]


    # minimize
    def print_guess(xk):
        print(f"Current Parameters: {xk}")

    # Example usage in minimize
    result = minimize(sumSq_all, 
                      initial_guess,
                      bounds=bounds,
                      callback=print_guess
                      )


    return result

def estimate_feedback():
    # set bounds
    a_f0 = a_f
    a_f_bounds = (-20e-4, -1.0e-4)

    a_b0 = -a_b
    a_b_bounds = (-5e-5, 5e-5)

    a_c0 = -a_c
    a_c_bounds = (0.5e-4, 20e-4)

    ins0 = 4e-3 
    ins_bounds = (1e-3, 20e-3)

    b0 = beta_t
    b_bounds = (0.5e-3,10e-3)

    initial_guess = [a_f0,a_b0,a_c0,ins0,b0]
    # print(initial_guess)
    bounds = [a_f_bounds,a_b_bounds,a_c_bounds,ins_bounds,b_bounds]

    # minimize
    result = minimize(sumSq_initial, initial_guess, bounds=bounds,method='Nelder-Mead')
    # result = 0
    return result

def main():
    # # %%
    result = estimate_ht()

    # # %%
    print(result.x)
    # Assuming 'result.x' contains the value you want to write to the file
    value_to_write = result.x

    # Specify the file path where you want to save the value
    file_path = "estimate_hA_all_2.txt"

    # Open the file in write mode and write the value to it
    with open(file_path, "w") as file:
        file.write(str(value_to_write))

    # The value has been written to the file
    print(f"Value has been written to {file_path}")

    # from_file = [0.02187095, 0.0113666,  0.00104331, 0.00415982, 0.00266767, 0.0093217, 0.00080961, 0.00687754, 0.19685921, 0.00049369, 0.0741425]
    # a_f0 = a_f
    # a_b0 = a_b
    # a_c0 = a_c
    # ins0 = 4e-3
    # b0 = beta_t
    # initial_guess = [a_f,a_b,a_c,ins0,beta_t]
    # res = relax_fib_oop(from_file)
    # simulation_output = [s[6]*P for s in res][i_insert[0]:(i_insert[-1]+1)]
    # # print(simulation_output[-1])
    # print(res[-1][6])
    # plt.plot(T_insert,simulation_output,label="JiTCDDE")
    # plt.plot(T_insert,interpolated_values,label="ORNL-1845")
    # plt.xlabel(r"$t$ (s)")
    # plt.ylabel("MW")
    # plt.title("Reactivity Insertions")
    # plt.legend()
    # plt.savefig("test_full_report.png")
    # plt.show()

if __name__ == '__main__':
    main()