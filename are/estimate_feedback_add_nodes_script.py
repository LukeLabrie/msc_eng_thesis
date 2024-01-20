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
# d = df[1][0]-P
# df[1] = [p-d for p in df[1]]

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
        sol_jit = relax_ht(params) # Implement this function to run your model
        
        # calculate error
        simulation_output = [s[46]*P for s in sol_jit][i_insert[0]:(i_insert[-1]+1)]
        error = sum((simulation_output - interpolated_values)**2)  # Sum of squared errors
        return error
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
        sol_jit = relax_feedback_insertion_beta(params)
        
        # calculate error
        simulation_output = [s[6]*P for s in sol_jit][i_insert[0]:(i_insert[-1]+1)]
        error = sum((simulation_output - interpolated_values)**2)  # Sum of squared errors
        return error
    except:
        return float('inf')
    
def sumSq_add_nodes(params):
    # try:
        # run
    sol_jit = relax_feedback_extra_nodes(params)
    
    # calculate error
    simulation_output = [s[6]*P for s in sol_jit][i_insert[0]:(i_insert[-1]+1)]
    error = sum((simulation_output - interpolated_values)**2)  # Sum of squared errors
    return error
    # except:
    #     return float('inf')
# %%

from scipy.optimize import minimize

def print_guess(xk):
    # Print the current parameters
    print(f"Current Parameters: {xk}")

    # Append the parameters to a text file
    with open('parameters_log.txt', 'a') as file:
        file.write(f"{xk}\n")


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


    initial_guess = [ft_c,tc_c,mc_c,ft_hx,ht_hx,ct_hx,th_hxch,ht_hxhw,
                     tw_hxhw,ht_hxhwc,tw_hxhwc]
    
    bounds = [ft_c_lims,tc_c_lims,mc_c_lims,ft_hx_lims,ht_hx_lims,ct_hx_lims,
              th_hxch_lims,ht_hxhw_lims,tw_hxhw_lims,ht_hxhwc_lims,tw_hxhwc_lims]

    # Example usage in minimize
    result = minimize(
        sumSq_all, 
        initial_guess, 
        bounds=bounds,
        callback=print_guess
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
    a_f_bounds = (-20e-5, 20e-5)

    a_b0 = -a_b
    a_b_bounds = (-5e-5, 5e-5)

    a_c0 = a_c
    a_c_bounds = (-20e-5, 20e-5)

    ins0 = inserted 
    ins_bounds = (1e-3, 20e-3)

    b0 = beta_t
    b_bounds = (1e-3,10e-3)

    initial_guess = [a_f0,a_b0,a_c0,ins0,b0]
    bounds = [a_f_bounds,a_b_bounds,a_c_bounds,ins_bounds,b_bounds]

    # minimize
    result = minimize(sumSq_initial, initial_guess, bounds=bounds,method='Nelder-Mead')

    return result

def estimate_feedback_add_nodes():
    # set bounds
    a_f0 = a_f
    a_f_bounds = (-20e-5, -0.5e-5)

    a_b0 = -a_b
    a_b_bounds = (-5e-5, 5e-5)

    a_c0 = a_c
    a_c_bounds = (-20e-5, 0.0)

    ins0 = inserted 
    ins_bounds = (0.5e-3, 20e-3)

    b0 = beta_t
    b_bounds = (1e-3,10e-3)

    initial_guess = [a_f0,a_b0,a_c0,ins0,b0]
    bounds = [a_f_bounds,a_b_bounds,a_c_bounds,ins_bounds,b_bounds]

    # minimize
    result = minimize(sumSq_add_nodes, 
                      initial_guess, 
                      bounds=bounds, 
                      method='Nelder-Mead',
                      callback=print_guess)

    return result


# %%
result = estimate_feedback_add_nodes()

# # # %%
print(result.x)
# # # Assuming 'result.x' contains the value you want to write to the file
value_to_write = result.x

# # # Specify the file path where you want to save the value
file_path = "output_feedback_add_nodes.txt"

# # # Open the file in write mode and write the value to it
with open(file_path, "w") as file:
    file.write(str(value_to_write))

# # # The value has been written to the file
print(f"Value has been written to {file_path}")

#from_file = [-1.99496857e-04, -4.99905872e-05,  1.74383022e-04,  6.82744637e-03, 3.52626522e-03] # no reinsertion
#from_file = [-2.36341143e-05, -5.00000000e-05,  2.01131763e-05,  1.00000000e-03, 2.97544233e-03] # reinsertion, forgot rho function
#from_file = [-4.12900591e-05,  4.33967181e-05,  1.85146114e-05,  1.76526592e-03, 2.22301472e-03]
#from_file = [-2.75109719e-05,  1.35671394e-05,  1.92674101e-05,  1.01544816e-03, 3.38047855e-03] #NM
# from_file = [-2.00000000e-04, 4.99999886e-05, 1.47679144e-04, 7.31105685e-03, 1.00000000e-03]
# # a_f0 = a_f
# # a_b0 = a_b
# # a_c0 = a_c
# # ins0 = 0.0
# # b0 = beta_t
# # initial_guess = [a_f,a_b,a_c,0.0,beta_t]
# res = relax_feedback_extra_nodes(from_file)
# simulation_output = [s[6]*P for s in res][i_insert[0]:(i_insert[-1]+1)]
# print(simulation_output-interpolated_values)

# print(simulation_output[-1])
# plt.plot(T_insert,simulation_output,label="JiTCDDE")
# plt.plot(T_insert,interpolated_values,label="ORNL-1845")
# plt.xlabel(r"$t$ (s)")
# plt.ylabel("MW")
# plt.title("Reactivity Insertions")
# plt.legend()
# plt.savefig("test_full_report.png")
# plt.show()