# %%
from relaxations import *
from parameters import *
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
from scipy.interpolate import CubicSpline

# %%
# read experimental data
# unpack data
df_inlet_reversed = pd.read_csv("./data/fuel_inlet_temp.csv",header=None)
df_outlet_reversed = pd.read_csv("./data/fuel_outlet_temp.csv",header=None)
df_inlet = df_inlet_reversed.iloc[::-1]
df_outlet = df_outlet_reversed.iloc[::-1]
df_inlet = df_inlet.reset_index(drop=True)
df_outlet = df_outlet.reset_index(drop=True)
df_inlet[1] = [F_to_K(t) for t in df_inlet[1]]
df_outlet[1] = [F_to_K(t) for t in df_outlet[1]]

# get indices for simulation data
t_before_data = (1110-df_inlet[0][0])*60
duration_data = (df_inlet.iloc[-1][0]-df_inlet[0][0])*60
t_end_data = df_inlet.iloc[-1][0]
t_before_sim = t_ins-t_before_data
T_insert = [t for t in T if (t > (t_before_sim)) and (t < (t_before_sim)+(duration_data))]
i_insert = [t[0] for t in enumerate(T) if (t[1] > (t_before_sim)) and (t[1] < (t_before_sim)+(duration_data))]
adj = (df_inlet[0][0])*60-T_insert[0]

# Set up interpolation
# Assuming df[0] is time and df[1] is the data you want to interpolate
spline_inlet = CubicSpline(df_inlet[0], df_inlet[1])  # Multiplying df[0] by 60 if it's in minutes
spline_outlet = CubicSpline(df_outlet[0], df_outlet[1])  # Multiplying df[0] by 60 if it's in minutes



# Use the spline to interpolate at the desired times
interpolated_inlet = spline_inlet(T[i_insert[0]:i_insert[-1]+1])
interpolated_outlet = spline_outlet(T[i_insert[0]:i_insert[-1]+1])


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
    
def sumSq_temp(params):
    # try:
    sol_jit = relax_temp(params) # Implement this function to run your model
    
    # calculate error
    simulation_output_inlet = [s[0] for s in sol_jit][i_insert[0]:(i_insert[-1]+1)]
    simulation_output_outlet = [s[1] for s in sol_jit][i_insert[0]:(i_insert[-1]+1)]
    error1 = sum((simulation_output_inlet - interpolated_inlet)**2)
    error2 = sum((simulation_output_outlet - interpolated_outlet)**2) 
    return error1 + error2
    # except:
    #      print("failed")
    #      return float('inf')
    

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

    def print_guess(xk):
        print(f"Current Parameters: {xk}")

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
    a_f_bounds = (-20e-5, -2e-5)

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


    # def write_guess(xk):
    #     f = 'inputs_all.txt'
    #     with open(f, 'a') as file:  # Open the file in append mode
    #         file.write(f"Current Parameters: {xk}\n")  # Write the parameters and a newline character

    # Example usage in minimize
    result = minimize(sumSq_all, 
                    initial_guess,
                    bounds=bounds,
                    method='Nelder-Mead'
                  )

    return result

def estimate_temp():
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

    scp_f0 = scp_f
    scp_f_lims = (scp_f/3,3*scp_f)

    scp_c0 = scp_c
    scp_c_lims = (scp_c/3,3*scp_c)

    initial_guess = [a_f0,a_b0,a_c0,ins0,b0,ft_c, tc_c,mc_c,ft_hx,ht_hx,ct_hx,
                     th_hxch,ht_hxhw,tw_hxhw,ht_hxhwc,tw_hxhwc,scp_f0,scp_c0]
    
    bounds = [a_f_bounds,a_b_bounds,a_c_bounds,ins_bounds,b_bounds, 
              ft_c_lims,tc_c_lims, mc_c_lims,ft_hx_lims,ht_hx_lims,ct_hx_lims,
              th_hxch_lims,ht_hxhw_lims,tw_hxhw_lims,ht_hxhwc_lims,tw_hxhwc_lims,
              scp_f_lims, scp_c_lims]
    
    result = minimize(sumSq_temp, 
                    initial_guess,
                    bounds=bounds,
                    method='Nelder-Mead'
                  )

    return result


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

scp_f0 = scp_f
scp_f_lims = (scp_f/3,3*scp_f)

scp_c0 = scp_c
scp_c_lims = (scp_c/3,3*scp_c)

initial_guess = [a_f0,a_b0,a_c0,ins0,b0,ft_c, tc_c,mc_c,ft_hx,ht_hx,ct_hx,
                    th_hxch,ht_hxhw,tw_hxhw,ht_hxhwc,tw_hxhwc,scp_f0,scp_c0]

result = sumSq_temp(initial_guess)
print(result)


# res = relax_all(res_NM)
# simulation_output = [s[6]*P for s in res][i_insert[0]:(i_insert[-1]+1)]
# plt.plot(T_insert,simulation_output,label="JiTCDDE")
# plt.plot(T_insert,interpolated_values,label="ORNL-1845")
# plt.xlabel(r"$t$ (s)")
# plt.ylabel("MW")
# plt.title("Reactivity Insertions")
# plt.legend()
# plt.savefig('test_all.png')
# plt.show()



# %%
# result = estimate_temp()

# # # %%
# print(result.x)
# # # Assuming 'result.x' contains the value you want to write to the file
# value_to_write = result.x

# # # Specify the file path where you want to save the value
# file_path = "output_estimate_temp.txt"

# # # Open the file in write mode and write the value to it
# with open(file_path, "w") as file:
#     file.write(str(value_to_write))

# # # The value has been written to the file
# print(f"Value has been written to {file_path}")

