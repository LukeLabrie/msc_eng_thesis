
from variables import *
from jitcdde import t, jitcdde_input, input
from heat_transfer import *
from parameters import *
from chspy import CubicHermiteSpline
from msrDynamics.objects import Node, System

def relax_feedback_insertion_beta(params):

    # Unpack parameters
    a_f, a_b, a_c, inserted, beta_t = params

    # CORE
    # fuel nodes
    dc_f1 = dT_bulkFlow(W_f, m_f_c/2, (hx_fh1_f2(t-tau_hx_c_f)+hx_fh2_f2(t-tau_hx_c_f))/2,c_f1()) + dT_internal(k_f1, P, mcp_f_c, n()) + F*dT_convective([c_t1()],c_f1(),[hA_ft_c/mcp_f_c]) - arbitrary_removal*c_f1()
    dc_f2 = dT_bulkFlow(W_f, m_f_c/2, c_f1(), c_f2()) +                                             dT_internal(k_f1, P, mcp_f_c, n()) + dT_convective([c_t1()],c_f2(),[hA_ft_c/mcp_f_c]) - arbitrary_removal*c_f2()

    # tubes
    dc_t1= F*dT_convective([c_f1(),c_f2(),c_c1(),c_c2()],c_t1(),[hA_ft_c/mcp_t_c,hA_ft_c/mcp_t_c,hA_tc_c/mcp_t_c,hA_tc_c/mcp_t_c]) - arbitrary_removal*c_t1()

    # coolant 
    dc_c1 = dT_bulkFlow(W_c,m_c_c/2,(hx_ch1_c2(t-tau_hx_c_c)+hx_ch2_c2(t-tau_hx_c_c))/2,c_c1()) + F*dT_convective([c_t1(),c_m1()],c_c1(),[hA_tc_c/mcp_c_c,hA_mc_c/mcp_c_c]) - arbitrary_removal*c_c1()
    dc_c2 = dT_bulkFlow(W_c,m_c_c/2,c_c1(),c_c2()) +                                              F*dT_convective([c_t1(),c_m1()],c_c2(),[hA_tc_c/mcp_c_c,hA_mc_c/mcp_c_c]) - arbitrary_removal*c_c2()

    # moderator 
    dc_m =  dT_internal(k_m,P,mcp_m_c,n()) + F*dT_convective([c_c1(),c_c2()],c_m1(),[hA_mc_c/mcp_m_c,hA_mc_c/mcp_m_c])

    # FUEL-HELIUM HX1
    dhx_fh1_f1 = dT_bulkFlow(W_f,m_f_hx,c_f2(t-tau_c_hx_f),hx_fh1_f1()) + F*dT_convective([hx_fh1_t1()],hx_fh1_f1(),[hA_ft_hx/mcp_f_hx]) 
    dhx_fh1_f2 = dT_bulkFlow(W_f,m_f_hx,hx_fh1_f1(),hx_fh1_f2()) +        F*dT_convective([hx_fh1_t1()],hx_fh1_f2(),[hA_ft_hx/mcp_f_hx]) 

    # Tubes 
    dhx_fh1_t1 = dT_convective([hx_fh1_f1(),hx_fh1_f2(),hx_fh1_h1(),hx_fh1_h2()],hx_fh1_t1(),[h/mcp_t_hx for h in [hA_ft_hx,hA_ft_hx,hA_ht_hx,hA_ht_hx]])

    # Helium
    dhx_fh1_h1 = dT_bulkFlow(W_h_fh,m_h_hxfh,hx_hwf2_h2(t-tau_h),hx_fh1_h1()) + F*dT_convective([hx_fh1_t1()],hx_fh1_h1(),[hA_ht_hx/mcp_h_hxfh])
    dhx_fh1_h2 = dT_bulkFlow(W_h_fh,m_h_hxfh,hx_fh1_h1(),hx_fh1_h2()) +         F*dT_convective([hx_fh1_t1()],hx_fh1_h2(),[hA_ht_hx/mcp_h_hxfh])

    # FUEL-HELIUM HX2 Fuel Nodes
    dhx_fh2_f1 = dT_bulkFlow(W_f,m_f_hx,c_f2(t-tau_c_hx_f), hx_fh2_f1()) + F*dT_convective([hx_fh2_t1()],hx_fh2_f1(),[hA_ft_hx/mcp_f_hx])
    dhx_fh2_f2 = dT_bulkFlow(W_f,m_f_hx,hx_fh2_f1(),hx_fh2_f2()) +         F*dT_convective([hx_fh2_t1()],hx_fh2_f2(),[hA_ft_hx/mcp_f_hx])

    # Tubes for FUEL-HELIUM HX2
    dhx_fh2_t1 = dT_convective([hx_fh2_f1(),hx_fh2_f2(),hx_fh2_h1(),hx_fh2_h2()],hx_fh2_t1(),[h/mcp_t_hx for h in [hA_ft_hx,hA_ft_hx,hA_ht_hx,hA_ht_hx]])

    # Helium for FUEL-HELIUM HX2
    dhx_fh2_h1 = dT_bulkFlow(W_h_fh,m_h_hxfh,hx_hwf1_h2(),hx_fh2_h1()) + F*dT_convective([hx_fh2_t1()],hx_fh2_h1(),[hA_ht_hx/mcp_h_hxfh])
    dhx_fh2_h2 = dT_bulkFlow(W_h_fh,m_h_hxfh,hx_fh2_h1(),hx_fh2_h2()) +  F*dT_convective([hx_fh2_t1()],hx_fh2_h2(),[hA_ht_hx/mcp_h_hxfh])

    # COOLANT-HELIUM HX1
    # Fuel Nodes
    dhx_ch1_c1 = dT_bulkFlow(W_c, m_c_hx, c_c2(t-tau_c_hx_f), hx_ch1_c1()) + F*dT_convective([hx_ch1_t1()], hx_ch1_c1(), [hA_ct_hx/mcp_c_hxch])
    dhx_ch1_c2 = dT_bulkFlow(W_c, m_c_hx, hx_ch1_c1(), hx_ch1_c2()) +        F*dT_convective([hx_ch1_t1()], hx_ch1_c2(), [hA_ct_hx/mcp_c_hxch])

    # Tubes
    dhx_ch1_t1 = F*dT_convective([hx_ch1_c1(),hx_ch1_c2(),hx_ch1_h1(),hx_ch1_h2()],hx_ch1_t1(),[h/mcp_t_hxch for h in [hA_ct_hx,hA_ct_hx,hA_th_hxch,hA_th_hxch]])

    # Helium
    dhx_ch1_h1 = dT_bulkFlow(W_h_ch, m_h_hxch/2, hx_hwc1_h2(t-tau_h), hx_ch1_h1()) + F*dT_convective([hx_ch1_t1()], hx_ch1_h1(), [hA_th_hxch/mcp_h_hxch])
    dhx_ch1_h2 = dT_bulkFlow(W_h_ch, m_h_hxch/2, hx_ch1_h1(), hx_ch1_h2()) +         F*dT_convective([hx_ch1_t1()], hx_ch1_h2(), [hA_th_hxch/mcp_h_hxch])

    # COOLANT-HELIUM HX2
    # Fuel Nodes
    dhx_ch2_c1 = dT_bulkFlow(W_c, m_c_hx, c_c2(t-tau_c_hx_f), hx_ch2_c1()) + F*dT_convective([hx_ch2_t1()], hx_ch2_c1(), [hA_ct_hx/mcp_c_hxch])
    dhx_ch2_c2 = dT_bulkFlow(W_c, m_c_hx, hx_ch2_c1(), hx_ch2_c2()) +        F*dT_convective([hx_ch2_t1()], hx_ch2_c2(), [hA_ct_hx/mcp_c_hxch])

    # Tubes
    dhx_ch2_t1 = dT_convective([hx_ch2_c1(),hx_ch2_c2(),hx_ch2_h1(),hx_ch2_h2()],hx_ch2_t1(),[h/mcp_t_hxch for h in [hA_ct_hx,hA_ct_hx,hA_th_hxch,hA_th_hxch]])

    # Helium
    dhx_ch2_h1 = dT_bulkFlow(W_h_ch, m_h_hxch/2, hx_hwc2_h2(t-tau_h), hx_ch2_h1()) + F*dT_convective([hx_ch2_t1()], hx_ch2_h1(), [hA_th_hxch/mcp_h_hxch])
    dhx_ch2_h2 = dT_bulkFlow(W_h_ch, m_h_hxch/2, hx_ch2_h1(), hx_ch2_h2()) +         F*dT_convective([hx_ch2_t1()], hx_ch2_h2(), [hA_th_hxch/mcp_h_hxch])

    # HELIUM-WATER HX1 (FUEL LOOP)
    # Helium
    T_hhwf_h1 = dT_bulkFlow(W_h_fh, m_h_hxhwf/2, hx_fh1_h2(), hx_hwf1_h1()) +  F*dT_convective([hx_hwf1_t1()], hx_hwf1_h1(), [hA_ht_hxhw/mcp_h_hxhw])
    T_hhwf_h2 = dT_bulkFlow(W_h_fh, m_h_hxhwf/2, hx_hwf1_h1(), hx_hwf1_h2()) + F*dT_convective([hx_hwf1_t1()], hx_hwf1_h2(), [hA_ht_hxhw/mcp_h_hxhw])

    # Tubes
    T_hhwf_t1 = F*dT_convective([hx_hwf1_h1(),hx_hwf1_h2(),hx_hwf1_w1(),hx_hwf1_w2()],hx_hwf1_t1(),[h/mcp_t_hxhw for h in [hA_ht_hxhw,hA_ht_hxhw,hA_tw_hxhw,hA_tw_hxhw]])

    # Water
    T_hhwf_w1 = dT_bulkFlow(W_hhwf_w, m_w_hxhwf/2, T0_hhwf_w1-hx_hwf1_w1(), hx_hwf1_w1(),dumped=True) + F*dT_convective([hx_hwf1_t1()], hx_hwf1_w1(), [hA_tw_hxhw/mcp_w])
    T_hhwf_w2 = dT_bulkFlow(W_hhwf_w, m_w_hxhwf/2, hx_hwf1_w1(), hx_hwf1_w2()) +                        F*dT_convective([hx_hwf1_t1()], hx_hwf1_w2(), [hA_tw_hxhw/mcp_w])

    # HELIUM-WATER HX2 (FUEL LOOP)
    # Helium
    T_hhwf2_h1 = dT_bulkFlow(W_h_fh, m_h_hxhwf/2, hx_fh2_h2(), hx_hwf2_h1()) +  F*dT_convective([hx_hwf2_t1()], hx_hwf2_h1(), [hA_ht_hxhw/mcp_h_hxhw])
    T_hhwf2_h2 = dT_bulkFlow(W_h_fh, m_h_hxhwf/2, hx_hwf2_h1(), hx_hwf2_h2()) + F*dT_convective([hx_hwf2_t1()], hx_hwf2_h2(), [hA_ht_hxhw/mcp_h_hxhw])

    # Tubes
    T_hhwf2_t1 = F*dT_convective([hx_hwf2_h1(),hx_hwf2_h2(),hx_hwf2_w1(),hx_hwf2_w2()],hx_hwf2_t1(),[h/mcp_t_hxhw for h in [hA_ht_hxhw,hA_ht_hxhw,hA_tw_hxhw,hA_tw_hxhw]])

    # Water
    T_hhwf2_w1 = dT_bulkFlow(W_hhwf_w, m_w_hxhwf/2, T0_hhwf_w1-hx_hwf2_w1(), hx_hwf2_w1(),dumped=True) + F*dT_convective([hx_hwf2_t1()], hx_hwf2_w1(), [hA_tw_hxhw/mcp_w])
    T_hhwf2_w2 = dT_bulkFlow(W_hhwf_w, m_w_hxhwf/2, hx_hwf2_w1(), hx_hwf2_w2()) +                        F*dT_convective([hx_hwf2_t1()], hx_hwf2_w2(), [hA_tw_hxhw/mcp_w])

    # HELIUM-WATER HX1 (COOLANT LOOP)
    # Helium
    T_hhwc_h1 = dT_bulkFlow(W_h_ch, m_h_hxhwc/2, hx_ch1_h2(), hx_hwc1_h1()) + F*dT_convective([hx_hwc1_t1()], hx_hwc1_h1(), [hA_ht_hxhwc/mcp_h_hxhwc])
    T_hhwc_h2 = dT_bulkFlow(W_h_ch, m_h_hxhwc/2, hx_hwc1_h1(), hx_hwc1_h2()) + F*dT_convective([hx_hwc1_t1()], hx_hwc1_h2(), [hA_ht_hxhwc/mcp_h_hxhwc])

    # Tubes
    T_hhwc_t1 = F*dT_convective([hx_hwc1_h1(),hx_hwc1_h2(),hx_hwc1_w1(),hx_hwc1_w2()],hx_hwc1_t1(),[h/mcp_t_hxhwc for h in [hA_ht_hxhwc,hA_ht_hxhwc,hA_tw_hxhwc,hA_tw_hxhwc]])

    # Water
    T_hhwc_w1 = dT_bulkFlow(W_hhwc_w, m_w_hxhwc/2, T0_hhwf_w1-hx_hwc1_w1(), hx_hwc1_w1(), dumped=True) + F*dT_convective([hx_hwc1_t1()], hx_hwc1_w1(), [hA_tw_hxhwc/mcp_w_hxhwc])
    T_hhwc_w2 = dT_bulkFlow(W_hhwc_w, m_w_hxhwc/2, hx_hwc1_w1(), hx_hwc1_w2()) + F*dT_convective([hx_hwc1_t1()], hx_hwc1_w2(), [hA_tw_hxhwc/mcp_w_hxhwc])


    # HELIUM-WATER HX2 (COOLANT LOOP)
    # Helium
    T_hhwc2_h1 = dT_bulkFlow(W_h_ch, m_h_hxhwc/2, hx_ch2_h2(), hx_hwc2_h1()) + F*dT_convective([hx_hwc2_t1()], hx_hwc2_h1(), [hA_ht_hxhwc/mcp_h_hxhwc])
    T_hhwc2_h2 = dT_bulkFlow(W_h_ch, m_h_hxhwc/2, hx_hwc2_h1(), hx_hwc2_h2()) + F*dT_convective([hx_hwc2_t1()], hx_hwc2_h2(), [hA_ht_hxhwc/mcp_h_hxhwc])

    # Tubes
    T_hhwc2_t1 = F*dT_convective([hx_hwc2_h1(),hx_hwc2_h2(),hx_hwc2_w1(),hx_hwc2_w2()],hx_hwc2_t1(),[h/mcp_t_hxhwc for h in [hA_ht_hxhwc,hA_ht_hxhwc,hA_tw_hxhwc,hA_tw_hxhwc]])

    # Water
    T_hhwc2_w1 = dT_bulkFlow(W_hhwc_w, m_w_hxhwc/2, T0_hhwf_w1-hx_hwc2_w1(), hx_hwc2_w1(), dumped=True) + F*dT_convective([hx_hwc2_t1()], hx_hwc2_w1(), [hA_tw_hxhwc/mcp_w_hxhwc])
    T_hhwc2_w2 = dT_bulkFlow(W_hhwc_w, m_w_hxhwc/2, hx_hwc2_w1(), hx_hwc2_w2()) + F*dT_convective([hx_hwc2_t1()], hx_hwc2_w2(), [hA_tw_hxhwc/mcp_w_hxhwc])

    # # reactivity insertion
    def rho_insert(t):
        if (t<t_ins):
            return 0.0
        elif (t<(t_ins+insert_duration)):
            return ((t-t_ins))*(inserted/insert_duration) # linear
        elif (t < t_wd):
            return inserted
        elif (t < t_wd+insert_duration):
            return inserted-((t-t_wd))*(inserted/insert_duration) # linear
        else:
            return 0.0

    rho_spline = CubicHermiteSpline(n=1)
    rho_spline.from_function(rho_insert, times_of_interest = T)
    rho_ext = input(0)

    # def rho_insert(t):
    #     insert_duration = 0.4/0.011
    #     if (t<t_ins):
    #         return 0.0
    #     elif (t<(t_ins+insert_duration)):
    #         return ((t-t_ins))*(inserted/insert_duration)
    #     else:
    #         return inserted

    # rho_spline = CubicHermiteSpline(n=1)
    # rho_spline.from_function(rho_insert, times_of_interest = T)
    # rho_ext = input(0)

    dn = ((rho()+rho_ext)-beta_t)*n()/Lam+lam[0]*C1()+lam[1]*C2()+lam[2]*C3()+lam[3]*C4()+lam[4]*C5()+lam[5]*C6()           # n (no source insertion): n()
    
    # dC_i/dt (precursor concentrations)
    dC1 = n()*beta_t*(beta/beta_t)[0]/Lam - lam[0]*C1() - C1()/tau_c + C1(t-tau_l)*np.exp(-lam[0]*tau_l)/tau_c                       # C1: y(27)
    dC2 = n()*beta_t*(beta/beta_t)[1]/Lam - lam[1]*C2() - C2()/tau_c + C2(t-tau_l)*np.exp(-lam[1]*tau_l)/tau_c                       # C2: y(28)
    dC3 = n()*beta_t*(beta/beta_t)[2]/Lam - lam[2]*C3() - C3()/tau_c + C3(t-tau_l)*np.exp(-lam[2]*tau_l)/tau_c                       # C3: y(29)
    dC4 = n()*beta_t*(beta/beta_t)[3]/Lam - lam[3]*C4() - C4()/tau_c + C4(t-tau_l)*np.exp(-lam[3]*tau_l)/tau_c                       # C4: y(30)
    dC5 = n()*beta_t*(beta/beta_t)[4]/Lam - lam[4]*C5() - C5()/tau_c + C5(t-tau_l)*np.exp(-lam[4]*tau_l)/tau_c                       # C5: y(31)
    dC6 = n()*beta_t*(beta/beta_t)[5]/Lam - lam[5]*C6() - C6()/tau_c + C6(t-tau_l)*np.exp(-lam[5]*tau_l)/tau_c                       # C6: y(32)

    # reactivity 
    drho = (a_f/2)*(dc_f1 + dc_f2)+(a_b)*(dc_m)+(a_c/2)*(dc_c1+dc_c2)           # rho: y(33)

    # instantiate jitcdde object
    DDE = jitcdde_input([dc_f1,dc_f2,dc_t1,dc_c1,dc_c2,dc_m,dn,dC1,dC2,dC3,dC4,dC5,dC6,drho,
                        dhx_fh1_f1,dhx_fh1_f2,dhx_fh1_t1,dhx_fh1_h1,dhx_fh1_h2,
                        dhx_fh2_f1,dhx_fh2_f2,dhx_fh2_t1,dhx_fh2_h1,dhx_fh2_h2,
                        dhx_ch1_c1,dhx_ch1_c2,dhx_ch1_t1,dhx_ch1_h1,dhx_ch1_h2,
                        dhx_ch2_c1,dhx_ch2_c2,dhx_ch2_t1,dhx_ch2_h1,dhx_ch2_h2,
                        T_hhwf_h1,T_hhwf_h2,T_hhwf_t1,T_hhwf_w1,T_hhwf_w2,
                        T_hhwf2_h1,T_hhwf2_h2,T_hhwf2_t1,T_hhwf2_w1,T_hhwf2_w2,
                        T_hhwc_h1,T_hhwc_h2,T_hhwc_t1,T_hhwc_w1,T_hhwc_w2,
                        T_hhwc2_h1,T_hhwc2_h2,T_hhwc2_t1,T_hhwc2_w1,T_hhwc2_w2],
                        rho_spline)

    # set initial conditions
    DDE.constant_past([T0_c_f1,T0_c_f2,T0_c_t1,T0_c_c1,T0_c_c2,T0_c_m+50,n_frac0,C0[0],C0[1],C0[2],C0[3],C0[4],C0[5],0.0,
                    T0_hfh_f1,T0_hfh_f2,T0_hfh_t1,T0_hfh_h1,T0_hfh_h2,
                    T0_hfh_f1,T0_hfh_f2,T0_hfh_t1,T0_hfh_h1,T0_hfh_h2,
                    T0_hch_c1,T0_hch_c2,T0_hch_t1,T0_hch_h1,T0_hch_h2,
                    T0_hch_c1,T0_hch_c2,T0_hch_t1,T0_hch_h1,T0_hch_h2,
                    T0_hhwf_h1,T0_hhwf_h2,T0_hhwf_t1,T0_hhwf_w1,T0_hhwf_w2,
                    T0_hhwf_h1,T0_hhwf_h2,T0_hhwf_t1,T0_hhwf_w1,T0_hhwf_w2,
                    T0_hhwc_h1,T0_hhwc_h2,T0_hhwc_t1,T0_hhwc_w1,T0_hhwc_w2,
                    T0_hhwc_h1,T0_hhwc_h2,T0_hhwc_t1,T0_hhwc_w1,T0_hhwc_w2
                        ])

    sol_jit = []
    for t_x in T:
        sol_jit.append(DDE.integrate(t_x))

    print([dc_f1,dc_f2,dc_t1,dc_c1,dc_c2,dc_m,dn,dC1,dC2,dC3,dC4,dC5,dC6,drho,
                        dhx_fh1_f1,dhx_fh1_f2,dhx_fh1_t1,dhx_fh1_h1,dhx_fh1_h2,
                        dhx_fh2_f1,dhx_fh2_f2,dhx_fh2_t1,dhx_fh2_h1,dhx_fh2_h2,
                        dhx_ch1_c1,dhx_ch1_c2,dhx_ch1_t1,dhx_ch1_h1,dhx_ch1_h2,
                        dhx_ch2_c1,dhx_ch2_c2,dhx_ch2_t1,dhx_ch2_h1,dhx_ch2_h2,
                        T_hhwf_h1,T_hhwf_h2,T_hhwf_t1,T_hhwf_w1,T_hhwf_w2,
                        T_hhwf2_h1,T_hhwf2_h2,T_hhwf2_t1,T_hhwf2_w1,T_hhwf2_w2,
                        T_hhwc_h1,T_hhwc_h2,T_hhwc_t1,T_hhwc_w1,T_hhwc_w2,
                        T_hhwc2_h1,T_hhwc2_h2,T_hhwc2_t1,T_hhwc2_w1,T_hhwc2_w2])
    return sol_jit

def relax_fib_correct(params):

    # Unpack parameters
    a_f, a_b, a_c, inserted, beta_t = params

    # CORE
    # fuel nodes
    dc_f1 = dT_bulkFlow(W_f, m_f_c/2, (hx_fh1_f2(t-tau_hx_c_f)+hx_fh2_f2(t-tau_hx_c_f))/2,c_f1()) + dT_internal(k_f1, P, mcp_f_c/2, n()) + dT_convective([c_t1()],c_f1(),[hA_ft_c/mcp_f_c]) 
    dc_f2 = dT_bulkFlow(W_f, m_f_c/2, c_f1(), c_f2()) +                                             dT_internal(k_f1, P, mcp_f_c/2, n()) + dT_convective([c_t1()],c_f2(),[hA_ft_c/mcp_f_c])

    # tubes
    dc_t1= dT_convective([c_f1(),c_f2(),c_c1(),c_c2()],c_t1(),[(hA_ft_c/2)/mcp_t_c,(hA_ft_c/2)/mcp_t_c,(hA_tc_c/2)/mcp_t_c,(hA_tc_c/2)/mcp_t_c])

    # coolant 
    dc_c1 = dT_bulkFlow(W_c,m_c_c/2,(hx_ch1_c2(t-tau_hx_c_c)+hx_ch2_c2(t-tau_hx_c_c))/2,c_c1()) + dT_convective([c_t1(),c_m1()],c_c1(),[hA_tc_c/mcp_c_c,hA_mc_c/mcp_c_c])
    dc_c2 = dT_bulkFlow(W_c,m_c_c/2,c_c1(),c_c2()) +                                              dT_convective([c_t1(),c_m1()],c_c2(),[hA_tc_c/mcp_c_c,hA_mc_c/mcp_c_c])

    # moderator 
    dc_m =  dT_internal(k_m,P,mcp_m_c,n()) + dT_convective([c_c1(),c_c2()],c_m1(),[(hA_mc_c/2)/mcp_m_c,(hA_mc_c/2)/mcp_m_c])

    # FUEL-HELIUM HX1
    dhx_fh1_f1 = dT_bulkFlow(W_f,m_f_hx,c_f2(t-tau_c_hx_f),hx_fh1_f1()) + dT_convective([hx_fh1_t1()],hx_fh1_f1(),[hA_ft_hx/mcp_f_hx]) 
    dhx_fh1_f2 = dT_bulkFlow(W_f,m_f_hx,hx_fh1_f1(),hx_fh1_f2()) +        dT_convective([hx_fh1_t1()],hx_fh1_f2(),[hA_ft_hx/mcp_f_hx]) 

    # Tubes 
    dhx_fh1_t1 = dT_convective([hx_fh1_f1(),hx_fh1_f2(),hx_fh1_h1(),hx_fh1_h2()],hx_fh1_t1(),[(h/2)/mcp_t_hx for h in [hA_ft_hx,hA_ft_hx,hA_ht_hx,hA_ht_hx]])

    # Helium
    dhx_fh1_h1 = dT_bulkFlow(W_h_fh,m_h_hxfh/2,hx_hwf2_h2(t-tau_h),hx_fh1_h1()) + dT_convective([hx_fh1_t1()],hx_fh1_h1(),[hA_ht_hx/mcp_h_hxfh])
    dhx_fh1_h2 = dT_bulkFlow(W_h_fh,m_h_hxfh/2,hx_fh1_h1(),hx_fh1_h2()) +         dT_convective([hx_fh1_t1()],hx_fh1_h2(),[hA_ht_hx/mcp_h_hxfh])

    # FUEL-HELIUM HX2 Fuel Nodes
    dhx_fh2_f1 = dT_bulkFlow(W_f,m_f_hx,c_f2(t-tau_c_hx_f), hx_fh2_f1()) + dT_convective([hx_fh2_t1()],hx_fh2_f1(),[hA_ft_hx/mcp_f_hx])
    dhx_fh2_f2 = dT_bulkFlow(W_f,m_f_hx,hx_fh2_f1(),hx_fh2_f2()) +         dT_convective([hx_fh2_t1()],hx_fh2_f2(),[hA_ft_hx/mcp_f_hx])

    # Tubes for FUEL-HELIUM HX2
    dhx_fh2_t1 = dT_convective([hx_fh2_f1(),hx_fh2_f2(),hx_fh2_h1(),hx_fh2_h2()],hx_fh2_t1(),[(h/2)/mcp_t_hx for h in [hA_ft_hx,hA_ft_hx,hA_ht_hx,hA_ht_hx]])

    # Helium for FUEL-HELIUM HX2
    dhx_fh2_h1 = dT_bulkFlow(W_h_fh,m_h_hxfh/2,hx_hwf1_h2(),hx_fh2_h1()) + dT_convective([hx_fh2_t1()],hx_fh2_h1(),[hA_ht_hx/mcp_h_hxfh])
    dhx_fh2_h2 = dT_bulkFlow(W_h_fh,m_h_hxfh/2,hx_fh2_h1(),hx_fh2_h2()) +  dT_convective([hx_fh2_t1()],hx_fh2_h2(),[hA_ht_hx/mcp_h_hxfh])

    # COOLANT-HELIUM HX1
    # Fuel Nodes
    dhx_ch1_c1 = dT_bulkFlow(W_c, m_c_hx, c_c2(t-tau_c_hx_f), hx_ch1_c1()) + dT_convective([hx_ch1_t1()], hx_ch1_c1(), [hA_ct_hx/mcp_c_hxch])
    dhx_ch1_c2 = dT_bulkFlow(W_c, m_c_hx, hx_ch1_c1(), hx_ch1_c2()) +        dT_convective([hx_ch1_t1()], hx_ch1_c2(), [hA_ct_hx/mcp_c_hxch])

    # Tubes
    dhx_ch1_t1 = F*dT_convective([hx_ch1_c1(),hx_ch1_c2(),hx_ch1_h1(),hx_ch1_h2()],hx_ch1_t1(),[(h/2)/mcp_t_hxch for h in [hA_ct_hx,hA_ct_hx,hA_th_hxch,hA_th_hxch]])

    # Helium
    dhx_ch1_h1 = dT_bulkFlow(W_h_ch, m_h_hxch/2, hx_hwc1_h2(t-tau_h), hx_ch1_h1()) + dT_convective([hx_ch1_t1()], hx_ch1_h1(), [hA_th_hxch/mcp_h_hxch])
    dhx_ch1_h2 = dT_bulkFlow(W_h_ch, m_h_hxch/2, hx_ch1_h1(), hx_ch1_h2()) +         dT_convective([hx_ch1_t1()], hx_ch1_h2(), [hA_th_hxch/mcp_h_hxch])

    # COOLANT-HELIUM HX2
    # Fuel Nodes
    dhx_ch2_c1 = dT_bulkFlow(W_c, m_c_hx, c_c2(t-tau_c_hx_f), hx_ch2_c1()) + dT_convective([hx_ch2_t1()], hx_ch2_c1(), [hA_ct_hx/mcp_c_hxch])
    dhx_ch2_c2 = dT_bulkFlow(W_c, m_c_hx, hx_ch2_c1(), hx_ch2_c2()) +        dT_convective([hx_ch2_t1()], hx_ch2_c2(), [hA_ct_hx/mcp_c_hxch])

    # Tubes
    dhx_ch2_t1 = dT_convective([hx_ch2_c1(),hx_ch2_c2(),hx_ch2_h1(),hx_ch2_h2()],hx_ch2_t1(),[(h/2)/mcp_t_hxch for h in [hA_ct_hx,hA_ct_hx,hA_th_hxch,hA_th_hxch]])

    # Helium
    dhx_ch2_h1 = dT_bulkFlow(W_h_ch, m_h_hxch/2, hx_hwc2_h2(t-tau_h), hx_ch2_h1()) + dT_convective([hx_ch2_t1()], hx_ch2_h1(), [hA_th_hxch/mcp_h_hxch])
    dhx_ch2_h2 = dT_bulkFlow(W_h_ch, m_h_hxch/2, hx_ch2_h1(), hx_ch2_h2()) +         dT_convective([hx_ch2_t1()], hx_ch2_h2(), [hA_th_hxch/mcp_h_hxch])

    # HELIUM-WATER HX1 (FUEL LOOP)
    # Helium
    T_hhwf_h1 = dT_bulkFlow(W_h_fh, m_h_hxhwf/2, hx_fh1_h2(), hx_hwf1_h1()) +  dT_convective([hx_hwf1_t1()], hx_hwf1_h1(), [hA_ht_hxhw/mcp_h_hxhw])
    T_hhwf_h2 = dT_bulkFlow(W_h_fh, m_h_hxhwf/2, hx_hwf1_h1(), hx_hwf1_h2()) + dT_convective([hx_hwf1_t1()], hx_hwf1_h2(), [hA_ht_hxhw/mcp_h_hxhw])

    # Tubes
    T_hhwf_t1 = F*dT_convective([hx_hwf1_h1(),hx_hwf1_h2(),hx_hwf1_w1(),hx_hwf1_w2()],hx_hwf1_t1(),[(h/2)/mcp_t_hxhw for h in [hA_ht_hxhw,hA_ht_hxhw,hA_tw_hxhw,hA_tw_hxhw]])

    # Water
    T_hhwf_w1 = dT_bulkFlow(W_hhwf_w, m_w_hxhwf/2, T0_hhwf_w1-hx_hwf1_w1(), hx_hwf1_w1(),dumped=True) + dT_convective([hx_hwf1_t1()], hx_hwf1_w1(), [hA_tw_hxhw/mcp_w])
    T_hhwf_w2 = dT_bulkFlow(W_hhwf_w, m_w_hxhwf/2, hx_hwf1_w1(), hx_hwf1_w2()) +                        dT_convective([hx_hwf1_t1()], hx_hwf1_w2(), [hA_tw_hxhw/mcp_w])

    # HELIUM-WATER HX2 (FUEL LOOP)
    # Helium
    T_hhwf2_h1 = dT_bulkFlow(W_h_fh, m_h_hxhwf/2, hx_fh2_h2(), hx_hwf2_h1()) +  dT_convective([hx_hwf2_t1()], hx_hwf2_h1(), [hA_ht_hxhw/mcp_h_hxhw])
    T_hhwf2_h2 = dT_bulkFlow(W_h_fh, m_h_hxhwf/2, hx_hwf2_h1(), hx_hwf2_h2()) + dT_convective([hx_hwf2_t1()], hx_hwf2_h2(), [hA_ht_hxhw/mcp_h_hxhw])

    # Tubes
    T_hhwf2_t1 = F*dT_convective([hx_hwf2_h1(),hx_hwf2_h2(),hx_hwf2_w1(),hx_hwf2_w2()],hx_hwf2_t1(),[(h/2)/mcp_t_hxhw for h in [hA_ht_hxhw,hA_ht_hxhw,hA_tw_hxhw,hA_tw_hxhw]])

    # Water
    T_hhwf2_w1 = dT_bulkFlow(W_hhwf_w, m_w_hxhwf/2, T0_hhwf_w1-hx_hwf2_w1(), hx_hwf2_w1(),dumped=True) + dT_convective([hx_hwf2_t1()], hx_hwf2_w1(), [hA_tw_hxhw/mcp_w])
    T_hhwf2_w2 = dT_bulkFlow(W_hhwf_w, m_w_hxhwf/2, hx_hwf2_w1(), hx_hwf2_w2()) +                        dT_convective([hx_hwf2_t1()], hx_hwf2_w2(), [hA_tw_hxhw/mcp_w])

    # HELIUM-WATER HX1 (COOLANT LOOP)
    # Helium
    T_hhwc_h1 = dT_bulkFlow(W_h_ch, m_h_hxhwc/2, hx_ch1_h2(), hx_hwc1_h1()) + dT_convective([hx_hwc1_t1()], hx_hwc1_h1(), [hA_ht_hxhwc/mcp_h_hxhwc])
    T_hhwc_h2 = dT_bulkFlow(W_h_ch, m_h_hxhwc/2, hx_hwc1_h1(), hx_hwc1_h2()) + dT_convective([hx_hwc1_t1()], hx_hwc1_h2(), [hA_ht_hxhwc/mcp_h_hxhwc])

    # Tubes
    T_hhwc_t1 = F*dT_convective([hx_hwc1_h1(),hx_hwc1_h2(),hx_hwc1_w1(),hx_hwc1_w2()],hx_hwc1_t1(),[(h/2)/mcp_t_hxhwc for h in [hA_ht_hxhwc,hA_ht_hxhwc,hA_tw_hxhwc,hA_tw_hxhwc]])

    # Water
    T_hhwc_w1 = dT_bulkFlow(W_hhwc_w, m_w_hxhwc/2, T0_hhwc_w1-hx_hwc1_w1(), hx_hwc1_w1(), dumped=True) + dT_convective([hx_hwc1_t1()], hx_hwc1_w1(), [hA_tw_hxhwc/mcp_w_hxhwc])
    T_hhwc_w2 = dT_bulkFlow(W_hhwc_w, m_w_hxhwc/2, hx_hwc1_w1(), hx_hwc1_w2()) + dT_convective([hx_hwc1_t1()], hx_hwc1_w2(), [hA_tw_hxhwc/mcp_w_hxhwc])


    # HELIUM-WATER HX2 (COOLANT LOOP)
    # Helium
    T_hhwc2_h1 = dT_bulkFlow(W_h_ch, m_h_hxhwc/2, hx_ch2_h2(), hx_hwc2_h1()) + dT_convective([hx_hwc2_t1()], hx_hwc2_h1(), [hA_ht_hxhwc/mcp_h_hxhwc])
    T_hhwc2_h2 = dT_bulkFlow(W_h_ch, m_h_hxhwc/2, hx_hwc2_h1(), hx_hwc2_h2()) + dT_convective([hx_hwc2_t1()], hx_hwc2_h2(), [hA_ht_hxhwc/mcp_h_hxhwc])

    # Tubes
    T_hhwc2_t1 = dT_convective([hx_hwc2_h1(),hx_hwc2_h2(),hx_hwc2_w1(),hx_hwc2_w2()],hx_hwc2_t1(),[(h/2)/mcp_t_hxhwc for h in [hA_ht_hxhwc,hA_ht_hxhwc,hA_tw_hxhwc,hA_tw_hxhwc]])

    # Water
    T_hhwc2_w1 = dT_bulkFlow(W_hhwc_w, m_w_hxhwc/2, T0_hhwc_w1-hx_hwc2_w1(), hx_hwc2_w1(), dumped=True) + dT_convective([hx_hwc2_t1()], hx_hwc2_w1(), [hA_tw_hxhwc/mcp_w_hxhwc])
    T_hhwc2_w2 = dT_bulkFlow(W_hhwc_w, m_w_hxhwc/2, hx_hwc2_w1(), hx_hwc2_w2()) + dT_convective([hx_hwc2_t1()], hx_hwc2_w2(), [hA_tw_hxhwc/mcp_w_hxhwc])

    # # reactivity insertion
    def rho_insert(t):
        if (t<t_ins):
            return 0.0
        elif (t<(t_ins+insert_duration)):
            return ((t-t_ins))*(inserted/insert_duration) # linear
        elif (t < t_wd):
            return inserted
        elif (t < t_wd+insert_duration):
            return inserted-((t-t_wd))*(inserted/insert_duration) # linear
        else:
            return 0.0

    rho_spline = CubicHermiteSpline(n=1)
    rho_spline.from_function(rho_insert, times_of_interest = T)
    rho_ext = input(0)

    # def rho_insert(t):
    #     insert_duration = 0.4/0.011
    #     if (t<t_ins):
    #         return 0.0
    #     elif (t<(t_ins+insert_duration)):
    #         return ((t-t_ins))*(inserted/insert_duration)
    #     else:
    #         return inserted

    # rho_spline = CubicHermiteSpline(n=1)
    # rho_spline.from_function(rho_insert, times_of_interest = T)
    # rho_ext = input(0)

    dn = ((rho()+rho_ext)-beta_t)*n()/Lam+lam[0]*C1()+lam[1]*C2()+lam[2]*C3()+lam[3]*C4()+lam[4]*C5()+lam[5]*C6()           # n (no source insertion): n()
    
    # dC_i/dt (precursor concentrations)
    dC1 = n()*beta_t*(beta/beta_t)[0]/Lam - lam[0]*C1() - C1()/tau_c + C1(t-tau_l)*np.exp(-lam[0]*tau_l)/tau_c                       # C1: y(27)
    dC2 = n()*beta_t*(beta/beta_t)[1]/Lam - lam[1]*C2() - C2()/tau_c + C2(t-tau_l)*np.exp(-lam[1]*tau_l)/tau_c                       # C2: y(28)
    dC3 = n()*beta_t*(beta/beta_t)[2]/Lam - lam[2]*C3() - C3()/tau_c + C3(t-tau_l)*np.exp(-lam[2]*tau_l)/tau_c                       # C3: y(29)
    dC4 = n()*beta_t*(beta/beta_t)[3]/Lam - lam[3]*C4() - C4()/tau_c + C4(t-tau_l)*np.exp(-lam[3]*tau_l)/tau_c                       # C4: y(30)
    dC5 = n()*beta_t*(beta/beta_t)[4]/Lam - lam[4]*C5() - C5()/tau_c + C5(t-tau_l)*np.exp(-lam[4]*tau_l)/tau_c                       # C5: y(31)
    dC6 = n()*beta_t*(beta/beta_t)[5]/Lam - lam[5]*C6() - C6()/tau_c + C6(t-tau_l)*np.exp(-lam[5]*tau_l)/tau_c                       # C6: y(32)

    # reactivity 
    drho = (a_f/2)*(dc_f1 + dc_f2)+(a_b)*(dc_m)+(a_c/2)*(dc_c1+dc_c2)           # rho: y(33)

    # instantiate jitcdde object
    DDE = jitcdde_input([dc_f1,dc_f2,dc_t1,dc_c1,dc_c2,dc_m,dn,dC1,dC2,dC3,dC4,dC5,dC6,drho,
                        dhx_fh1_f1,dhx_fh1_f2,dhx_fh1_t1,dhx_fh1_h1,dhx_fh1_h2,
                        dhx_fh2_f1,dhx_fh2_f2,dhx_fh2_t1,dhx_fh2_h1,dhx_fh2_h2,
                        dhx_ch1_c1,dhx_ch1_c2,dhx_ch1_t1,dhx_ch1_h1,dhx_ch1_h2,
                        dhx_ch2_c1,dhx_ch2_c2,dhx_ch2_t1,dhx_ch2_h1,dhx_ch2_h2,
                        T_hhwf_h1,T_hhwf_h2,T_hhwf_t1,T_hhwf_w1,T_hhwf_w2,
                        T_hhwf2_h1,T_hhwf2_h2,T_hhwf2_t1,T_hhwf2_w1,T_hhwf2_w2,
                        T_hhwc_h1,T_hhwc_h2,T_hhwc_t1,T_hhwc_w1,T_hhwc_w2,
                        T_hhwc2_h1,T_hhwc2_h2,T_hhwc2_t1,T_hhwc2_w1,T_hhwc2_w2],
                        rho_spline)

    # set initial conditions
    DDE.constant_past([T0_c_f1,T0_c_f2,T0_c_t1,T0_c_c1,T0_c_c2,T0_c_m+50,n_frac0,C0[0],C0[1],C0[2],C0[3],C0[4],C0[5],0.0,
                    T0_hfh_f1,T0_hfh_f2,T0_hfh_t1,T0_hfh_h1,T0_hfh_h2,
                    T0_hfh_f1,T0_hfh_f2,T0_hfh_t1,T0_hfh_h1,T0_hfh_h2,
                    T0_hch_c1,T0_hch_c2,T0_hch_t1,T0_hch_h1,T0_hch_h2,
                    T0_hch_c1,T0_hch_c2,T0_hch_t1,T0_hch_h1,T0_hch_h2,
                    T0_hhwf_h1,T0_hhwf_h2,T0_hhwf_t1,T0_hhwf_w1,T0_hhwf_w2,
                    T0_hhwf_h1,T0_hhwf_h2,T0_hhwf_t1,T0_hhwf_w1,T0_hhwf_w2,
                    T0_hhwc_h1,T0_hhwc_h2,T0_hhwc_t1,T0_hhwc_w1,T0_hhwc_w2,
                    T0_hhwc_h1,T0_hhwc_h2,T0_hhwc_t1,T0_hhwc_w1,T0_hhwc_w2
                        ])

    sol_jit = []
    for t_x in T:
        sol_jit.append(DDE.integrate(t_x))

    return sol_jit

def relax_fib_oop(params):

    # Unpack parameters
    a_f, a_b, a_c, inserted, beta_t = params

    # ARE system        
    ARE = System()

    # CORE NODES
    c_f1 = Node(m = m_f_c/2, scp = scp_f, W = W_f, y0 = T0_c_f1)
    c_f2 = Node(m = m_f_c/2, scp = scp_f, W = W_f, y0 = T0_c_f2)
    c_t1 = Node(m = m_t, scp = scp_t, y0 = T0_c_t1)
    c_c1 = Node(m = m_c_c/2, scp = scp_c, W = W_c, y0 = T0_c_c1)
    c_c2 = Node(m = m_c_c/2, scp = scp_c, W = W_c, y0 = T0_c_c2) 
    c_m1 = Node(m = m_m_c, scp = scp_m, y0 = T0_c_m+50)
    n = Node(y0 = n_frac0)
    C1 = Node(y0 = C0[0])
    C2 = Node(y0 = C0[1])
    C3 = Node(y0 = C0[2])
    C4 = Node(y0 = C0[3])
    C5 = Node(y0 = C0[4])
    C6 = Node(y0 = C0[5])
    rho = Node(y0 = 0.0)

    # add reactivity input
    def rho_insert(t):
        if (t<t_ins):
            return 0.0
        elif (t<(t_ins+insert_duration)):
            return ((t-t_ins))*(inserted/insert_duration) # linear
        elif (t < t_wd):
            return inserted
        elif (t < t_wd+insert_duration):
            return inserted-((t-t_wd))*(inserted/insert_duration) # linear
        else:
            return 0.0

    rho_ext = ARE.add_input(rho_insert, T)

    # FUEL-HELIUM HX1
    hx_fh1_f1 = Node(m = m_f_hx/2, scp = scp_f, W = W_f/2, y0 = T0_hfh_f1)
    hx_fh1_f2 = Node(m = m_f_hx/2, scp = scp_f, W = W_f/2, y0 = T0_hfh_f2)
    hx_fh1_t1 = Node(m = m_t_hxfh, scp = scp_t, y0 = T0_hfh_t1)
    hx_fh1_h1 = Node(m = m_h_hxfh/2, scp = scp_h, W = W_h_fh, y0 = T0_hfh_h1)
    hx_fh1_h2 = Node(m = m_h_hxfh/2, scp = scp_h, W = W_h_fh, y0 = T0_hfh_h2)

    # FUEL-HELIUM HX2
    hx_fh2_f1 = Node(m = m_f_hx/2, scp = scp_f, W = W_f/2, y0 = T0_hfh_f1)
    hx_fh2_f2 = Node(m = m_f_hx/2, scp = scp_f, W = W_f/2, y0 = T0_hfh_f2)
    hx_fh2_t1 = Node(m = m_t_hxfh, scp = scp_t, y0 = T0_hfh_t1)
    hx_fh2_h1 = Node(m = m_h_hxfh/2, scp = scp_h, W = W_h_fh, y0 = T0_hfh_h1)
    hx_fh2_h2 = Node(m = m_h_hxfh/2, scp = scp_h, W = W_h_fh, y0 = T0_hfh_h2)

    # COOLANT-HELIUM HX1
    hx_ch1_c1 = Node(m = m_c_hx/2, scp = scp_c, W = W_c/2, y0 = T0_hch_c1)
    hx_ch1_c2 = Node(m = m_c_hx/2, scp = scp_c, W = W_c/2, y0 = T0_hch_c2)
    hx_ch1_t1 = Node(m = m_t_hxch, scp = scp_t, y0 = T0_hch_t1)
    hx_ch1_h1 = Node(m = m_h_hxch/2, scp = scp_h, W = W_h_ch, y0 = T0_hch_h1)
    hx_ch1_h2 = Node(m = m_h_hxch/2, scp = scp_h, W = W_h_ch, y0 = T0_hch_h2)

    # COOLANT-HELIUM HX2
    hx_ch2_c1 = Node(m = m_c_hx/2, scp = scp_c, W = W_c/2, y0 = T0_hch_c1)
    hx_ch2_c2 = Node(m = m_c_hx/2, scp = scp_c, W = W_c/2, y0 = T0_hch_c2)
    hx_ch2_t1 = Node(m = m_t_hxch, scp = scp_t, y0 = T0_hfh_t1)
    hx_ch2_h1 = Node(m = m_h_hxch/2, scp = scp_h, W = W_h_ch, y0 = T0_hch_h1)
    hx_ch2_h2 = Node(m = m_h_hxch/2, scp = scp_h, W = W_h_ch, y0 = T0_hch_h2)

    # HELIUM-WATER HX1 (FUEL LOOP)
    hx_hwf1_h1 = Node(m = m_h_hxhwf/2, scp = scp_h, W = W_h_fh, y0 = T0_hhwf_h1)
    hx_hwf1_h2 = Node(m = m_h_hxhwf/2, scp = scp_h, W = W_h_fh, y0 = T0_hhwf_h2)
    hx_hwf1_t1 = Node(m = m_t_hxhwf, scp = scp_t, y0 = T0_hhwf_t1)
    hx_hwf1_w1 = Node(m = m_w_hxhwf/2, scp = scp_w, W = W_hhwf_w, y0 = T0_hhwf_w1)
    hx_hwf1_w2 = Node(m = m_w_hxhwf/2, scp = scp_w, W = W_hhwf_w, y0 = T0_hhwf_w2)

    # HELIUM-WATER HX2 (FUEL LOOP)
    hx_hwf2_h1 = Node(m = m_h_hxhwf/2, scp = scp_h, W = W_h_fh, y0 = T0_hhwf_h1)
    hx_hwf2_h2 = Node(m = m_h_hxhwf/2, scp = scp_h, W = W_h_fh, y0 = T0_hhwf_h2)
    hx_hwf2_t1 = Node(m = m_t_hxhwf, scp = scp_t, y0 = T0_hhwf_t1)
    hx_hwf2_w1 = Node(m = m_w_hxhwf/2, scp = scp_w, W = W_hhwf_w, y0 = T0_hhwf_w1)
    hx_hwf2_w2 = Node(m = m_w_hxhwf/2, scp = scp_w, W = W_hhwf_w, y0 = T0_hhwf_w2)

    # HELIUM-WATER HX1 (COOLANT LOOP)
    hx_hwc1_h1 = Node(m = m_h_hxhwc/2, scp = scp_h, W = W_h_ch, y0 = T0_hhwc_h1)
    hx_hwc1_h2 = Node(m = m_h_hxhwc/2, scp = scp_h, W = W_h_ch, y0 = T0_hhwc_h2)
    hx_hwc1_t1 = Node(m = m_t_hxhwc, scp = scp_t, y0 = T0_hhwf_t1)
    hx_hwc1_w1 = Node(m = m_w_hxhwc/2, scp = scp_w, W = W_hhwc_w, y0 = T0_hhwc_w1)
    hx_hwc1_w2 = Node(m = m_w_hxhwc/2, scp = scp_w, W = W_hhwc_w, y0 = T0_hhwc_w2)

    # HELIUM-WATER HX2 (COOLANT LOOP)
    hx_hwc2_h1 = Node(m = m_h_hxhwc/2, scp = scp_h, W = W_h_ch, y0 = T0_hhwc_h1)
    hx_hwc2_h2 = Node(m = m_h_hxhwc/2, scp = scp_h, W = W_h_ch, y0 = T0_hhwc_h2)
    hx_hwc2_t1 = Node(m = m_t_hxhwc, scp = scp_t, y0 = T0_hhwf_t1)
    hx_hwc2_w1 = Node(m = m_w_hxhwc/2, scp = scp_w, W = W_hhwc_w, y0 = T0_hhwc_w1)
    hx_hwc2_w2 = Node(m = m_w_hxhwc/2, scp = scp_w, W = W_hhwc_w, y0 = T0_hhwc_w2)



    ARE.add_nodes([c_f1,c_f2,c_t1,c_c1,c_c2,c_m1,n,C1,C2,C3,C4,C5,C6,rho,
                hx_fh1_f1,hx_fh1_f2,hx_fh1_t1,hx_fh1_h1,hx_fh1_h2,
                hx_fh2_f1,hx_fh2_f2,hx_fh2_t1,hx_fh2_h1,hx_fh2_h2,
                hx_ch1_c1,hx_ch1_c2,hx_ch1_t1,hx_ch1_h1,hx_ch1_h2,
                hx_ch2_c1,hx_ch2_c2,hx_ch2_t1,hx_ch2_h1,hx_ch2_h2,
                hx_hwf1_h1,hx_hwf1_h2,hx_hwf1_t1,hx_hwf1_w1,hx_hwf1_w2,
                hx_hwf2_h1,hx_hwf2_h2,hx_hwf2_t1,hx_hwf2_w1,hx_hwf2_w2,
                hx_hwc1_h1,hx_hwc1_h2,hx_hwc1_t1,hx_hwc1_w1,hx_hwc1_w2,
                hx_hwc2_h1,hx_hwc2_h2,hx_hwc2_t1,hx_hwc2_w1,hx_hwc2_w2,
                ])
    
    # CORE
    c_f1.set_dTdt_bulkFlow(source = (hx_fh1_f2.y(t-tau_hx_c_f)+hx_fh2_f2.y(t-tau_hx_c_f))/2) 
    c_f1.set_dTdt_internal(source = n.y(), k = k_f1*P)
    c_f1.set_dTdt_convective(source = [c_t1.y()], hA = [hA_ft_c/2])

    c_f2.set_dTdt_bulkFlow(source = c_f1.y()) 
    c_f2.set_dTdt_internal(source = n.y(), k = k_f2*P)
    c_f2.set_dTdt_convective(source = [c_t1.y()], hA = [hA_ft_c/2])

    c_t1.set_dTdt_convective(source = [c_f1.y(), c_f2.y(), c_c1.y(), c_c2.y()], 
                            hA = [hA_ft_c/2,hA_ft_c/2,hA_tc_c/2,hA_tc_c/2])

    c_c1.set_dTdt_bulkFlow(source = (hx_ch1_c2.y(t-tau_c_hx_f)+hx_ch2_c2.y(t-tau_c_hx_f))/2)
    c_c1.set_dTdt_convective(source = [c_t1.y(),c_m1.y()], hA = [hA_tc_c/2,hA_mc_c/2])

    c_c2.set_dTdt_bulkFlow(source = c_c1.y())
    c_c2.set_dTdt_convective(source = [c_t1.y(),c_m1.y()], hA = [hA_tc_c/2,hA_mc_c/2])

    c_m1.set_dTdt_internal(source = n.y(), k = k_m*P)
    c_m1.set_dTdt_convective(source = [c_c1.y(), c_c2.y()], hA = [hA_mc_c/2]*2)

    n.set_dndt(rho.y()+rho_ext, beta_t, Lam, lam, [C1.y(), C2.y(), C3.y(), C4.y(), C5.y(), C6.y()])
    C1.set_dcdt(n.y(), beta[0], Lam, lam[0], tau_c, tau_l)
    C2.set_dcdt(n.y(), beta[1], Lam, lam[1], tau_c, tau_l)
    C3.set_dcdt(n.y(), beta[2], Lam, lam[2], tau_c, tau_l)
    C4.set_dcdt(n.y(), beta[3], Lam, lam[3], tau_c, tau_l)
    C5.set_dcdt(n.y(), beta[4], Lam, lam[4], tau_c, tau_l)
    C6.set_dcdt(n.y(), beta[5], Lam, lam[5], tau_c, tau_l)
    rho.set_drdt([c_f1.dydt(),c_f2.dydt(),c_m1.dydt(),c_c1.dydt(),c_c2.dydt()],[a_f/2,a_f/2,a_b,a_c/2,a_c/2])

    # FUEL-HELIUM HX1
    hx_fh1_f1.set_dTdt_bulkFlow(source = c_f2.y(t-tau_c_hx_f))
    hx_fh1_f1.set_dTdt_convective(source = [hx_fh1_t1.y()], hA = [hA_ft_hx/2])

    hx_fh1_f2.set_dTdt_bulkFlow(source = hx_fh1_f1.y())
    hx_fh1_f2.set_dTdt_convective(source = [hx_fh1_t1.y()], hA = [hA_ft_hx/2])

    hx_fh1_t1.set_dTdt_convective(source = [hx_fh1_f1.y(),hx_fh1_f2.y(),hx_fh1_h1.y(),hx_fh1_h2.y()],
                                hA = [hA_ft_hx/2,hA_ft_hx/2,hA_ht_hx/2,hA_ht_hx/2])

    hx_fh1_h1.set_dTdt_bulkFlow(source = hx_hwf2_h2.y(t-tau_h))
    hx_fh1_h1.set_dTdt_convective(source = [hx_fh1_t1.y()], hA = [hA_ht_hx/2])

    hx_fh1_h2.set_dTdt_bulkFlow(source = hx_fh1_h1.y())
    hx_fh1_h2.set_dTdt_convective(source = [hx_fh1_t1.y()], hA = [hA_ht_hx/2])

    # FUEL-HELIUM HX2
    hx_fh2_f1.set_dTdt_bulkFlow(source = c_f2.y(t-tau_c_hx_f))
    hx_fh2_f1.set_dTdt_convective(source = [hx_fh2_t1.y()], hA = [hA_ft_hx/2])

    hx_fh2_f2.set_dTdt_bulkFlow(source = hx_fh2_f1.y())
    hx_fh2_f2.set_dTdt_convective(source = [hx_fh2_t1.y()], hA = [hA_ft_hx/2])

    hx_fh2_t1.set_dTdt_convective(source = [hx_fh2_f1.y(),hx_fh2_f2.y(),hx_fh2_h1.y(),hx_fh2_h2.y()],
                                hA = [hA_ft_hx/2,hA_ft_hx/2,hA_ht_hx/2,hA_ht_hx/2])

    hx_fh2_h1.set_dTdt_bulkFlow(source = hx_hwf1_h2.y())
    hx_fh2_h1.set_dTdt_convective(source = [hx_fh2_t1.y()], hA = [hA_ht_hx/2])

    hx_fh2_h2.set_dTdt_bulkFlow(source = hx_fh2_h1.y())
    hx_fh2_h2.set_dTdt_convective(source = [hx_fh2_t1.y()], hA = [hA_ht_hx/2])

    # COOLANT-HELIUM HX1
    hx_ch1_c1.set_dTdt_bulkFlow(source = c_c2.y(t-tau_c_hx_f))
    hx_ch1_c1.set_dTdt_convective(source = [hx_ch1_t1.y()], hA = [hA_ct_hx/2])

    hx_ch1_c2.set_dTdt_bulkFlow(source = hx_ch1_c1.y())
    hx_ch1_c2.set_dTdt_convective(source = [hx_ch1_t1.y()], hA = [hA_ct_hx/2])

    hx_ch1_t1.set_dTdt_convective(source = [hx_ch1_c1.y(),hx_ch1_c2.y(),hx_ch1_h1.y(),hx_ch1_h2.y()],
                                hA = [hA_ct_hx/2,hA_ct_hx/2,hA_th_hxch/2,hA_th_hxch/2])

    hx_ch1_h1.set_dTdt_bulkFlow(source = hx_hwc1_h2.y(t-tau_h))
    hx_ch1_h1.set_dTdt_convective(source = [hx_ch1_t1.y()], hA = [hA_th_hxch/2])

    hx_ch1_h2.set_dTdt_bulkFlow(source = hx_ch1_h1.y())
    hx_ch1_h2.set_dTdt_convective(source = [hx_ch1_t1.y()], hA = [hA_th_hxch/2])

    # COOLANT-HELIUM HX2
    hx_ch2_c1.set_dTdt_bulkFlow(source = c_c2.y(t-tau_c_hx_f))
    hx_ch2_c1.set_dTdt_convective(source = [hx_ch2_t1.y()], hA = [hA_ct_hx/2])

    hx_ch2_c2.set_dTdt_bulkFlow(source = hx_ch2_c1.y())
    hx_ch2_c2.set_dTdt_convective(source = [hx_ch2_t1.y()], hA = [hA_ct_hx/2])

    hx_ch2_t1.set_dTdt_convective(source = [hx_ch2_c1.y(),hx_ch2_c2.y(),hx_ch2_h1.y(),hx_ch2_h2.y()],
                                hA = [hA_ct_hx/2,hA_ct_hx/2,hA_th_hxch/2,hA_th_hxch/2])

    hx_ch2_h1.set_dTdt_bulkFlow(source = hx_hwc2_h2.y(t-tau_h))
    hx_ch2_h1.set_dTdt_convective(source = [hx_ch2_t1.y()], hA = [hA_th_hxch/2])

    hx_ch2_h2.set_dTdt_bulkFlow(source = hx_ch2_h1.y())
    hx_ch2_h2.set_dTdt_convective(source = [hx_ch2_t1.y()], hA = [hA_th_hxch/2])

    # HELIUM-WATER HX1 (FUEL LOOP)
    hx_hwf1_h1.set_dTdt_bulkFlow(source = hx_fh1_h2.y())
    hx_hwf1_h1.set_dTdt_convective(source = [hx_hwf1_t1.y()], hA = [hA_ht_hxhw/2])

    hx_hwf1_h2.set_dTdt_bulkFlow(source = hx_hwf1_h1.y())
    hx_hwf1_h2.set_dTdt_convective(source = [hx_hwf1_t1.y()], hA = [hA_ht_hxhw/2])

    hx_hwf1_t1.set_dTdt_convective(source = [hx_hwf1_h1.y(),hx_hwf1_h2.y(),hx_hwf1_w1.y(),hx_hwf1_w2.y()],
                                hA = [hA_ht_hxhw/2,hA_ht_hxhw/2,hA_tw_hxhw/2,hA_tw_hxhw/2])

    hx_hwf1_w1.set_dTdt_bulkFlow(source = T0_hhwf_w1)
    hx_hwf1_w1.set_dTdt_convective(source = [hx_hwf1_t1.y()], hA = [hA_tw_hxhw/2])

    hx_hwf1_w2.set_dTdt_bulkFlow(source = hx_hwf1_w1.y())
    hx_hwf1_w2.set_dTdt_convective(source = [hx_hwf1_t1.y()], hA = [hA_tw_hxhw/2])

    # HELIUM-WATER HX2 (FUEL LOOP)
    hx_hwf2_h1.set_dTdt_bulkFlow(source = hx_fh2_h2.y())
    hx_hwf2_h1.set_dTdt_convective(source = [hx_hwf2_t1.y()], hA = [hA_ht_hxhw/2])

    hx_hwf2_h2.set_dTdt_bulkFlow(source = hx_hwf2_h1.y())
    hx_hwf2_h2.set_dTdt_convective(source = [hx_hwf2_t1.y()], hA = [hA_ht_hxhw/2])

    hx_hwf2_t1.set_dTdt_convective(source = [hx_hwf2_h1.y(),hx_hwf2_h2.y(),hx_hwf2_w1.y(),hx_hwf2_w2.y()],
                                hA = [hA_ht_hxhw/2,hA_ht_hxhw/2,hA_tw_hxhw/2,hA_tw_hxhw/2])

    hx_hwf2_w1.set_dTdt_bulkFlow(source = T0_hhwf_w1)
    hx_hwf2_w1.set_dTdt_convective(source = [hx_hwf2_t1.y()], hA = [hA_tw_hxhw/2])

    hx_hwf2_w2.set_dTdt_bulkFlow(source = hx_hwf2_w1.y())
    hx_hwf2_w2.set_dTdt_convective(source = [hx_hwf2_t1.y()], hA = [hA_tw_hxhw/2])

    # HELIUM-WATER HX1 (COOLANT LOOP)
    hx_hwc1_h1.set_dTdt_bulkFlow(source = hx_ch1_h2.y())
    hx_hwc1_h1.set_dTdt_convective(source = [hx_hwc1_t1.y()], hA = [hA_ht_hxhwc/2])

    hx_hwc1_h2.set_dTdt_bulkFlow(source = hx_hwc1_h1.y())
    hx_hwc1_h2.set_dTdt_convective(source = [hx_hwc1_t1.y()], hA = [hA_ht_hxhwc/2])

    hx_hwc1_t1.set_dTdt_convective(source = [hx_hwc1_h1.y(),hx_hwc1_h2.y(),hx_hwc1_w1.y(),hx_hwc1_w2.y()],
                                hA = [hA_ht_hxhwc/2,hA_ht_hxhwc/2,hA_tw_hxhwc/2,hA_tw_hxhwc/2])

    hx_hwc1_w1.set_dTdt_bulkFlow(source = T0_hhwc_w1)
    hx_hwc1_w1.set_dTdt_convective(source = [hx_hwc1_t1.y()], hA = [hA_tw_hxhwc/2])

    hx_hwc1_w2.set_dTdt_bulkFlow(source = hx_hwc1_w1.y())
    hx_hwc1_w2.set_dTdt_convective(source = [hx_hwc1_t1.y()], hA = [hA_tw_hxhwc/2])

    # HELIUM-WATER HX2 (COOLANT LOOP)
    hx_hwc2_h1.set_dTdt_bulkFlow(source = hx_ch2_h2.y())
    hx_hwc2_h1.set_dTdt_convective(source = [hx_hwc2_t1.y()], hA = [hA_ht_hxhwc/2])

    hx_hwc2_h2.set_dTdt_bulkFlow(source = hx_hwc2_h1.y())
    hx_hwc2_h2.set_dTdt_convective(source = [hx_hwc2_t1.y()], hA = [hA_ht_hxhwc/2])

    hx_hwc2_t1.set_dTdt_convective(source = [hx_hwc2_h1.y(),hx_hwc2_h2.y(),hx_hwc2_w1.y(),hx_hwc2_w2.y()],
                                hA = [hA_ht_hxhwc/2,hA_ht_hxhwc/2,hA_tw_hxhwc/2,hA_tw_hxhwc/2])

    hx_hwc2_w1.set_dTdt_bulkFlow(source = T0_hhwc_w1)
    hx_hwc2_w1.set_dTdt_convective(source = [hx_hwc2_t1.y()], hA = [hA_tw_hxhwc/2])

    hx_hwc2_w2.set_dTdt_bulkFlow(source = hx_hwc2_w1.y())
    hx_hwc2_w2.set_dTdt_convective(source = [hx_hwc2_t1.y()], hA = [hA_tw_hxhwc/2])

    ARE.finalize()
    DDE = ARE.integrator
    sol_jit = []
    for t_x in T:
        sol_jit.append(DDE.integrate(t_x))

    # force memory clearing 
    for node in ARE.nodes:
        node = None
    DDE = None 
    ARE = None 
    
    return sol_jit

def relax_delays(params):

    # Unpack parameters
    tau_hx_c_f, tau_hx_c_c, tau_c_hx_f, tau_h_loop_f, tau_h_btw_f, tau_c_hx_c,\
    tau_h_loop_c, tau_h_btw_c  = params

    tau_l = tau_c_hx_f + tau_hx_c_f

    # CORE
    # fuel nodes
    dc_f1 = dT_bulkFlow(W_f, m_f_c/2, (hx_fh1_f2(t-tau_hx_c_f)+hx_fh2_f2(t-tau_hx_c_f))/2,c_f1()) + dT_internal(k_f1, P, mcp_f_c, n()) + F*dT_convective([c_t1()],c_f1(),[hA_ft_c/mcp_f_c]) - arbitrary_removal*c_f1()
    dc_f2 = dT_bulkFlow(W_f, m_f_c/2, c_f1(), c_f2()) +                                             dT_internal(k_f1, P, mcp_f_c, n()) + dT_convective([c_t1()],c_f2(),[hA_ft_c/mcp_f_c]) - arbitrary_removal*c_f2()

    # tubes
    dc_t1= F*dT_convective([c_f1(),c_f2(),c_c1(),c_c2()],c_t1(),[hA_ft_c/mcp_t_c,hA_ft_c/mcp_t_c,hA_tc_c/mcp_t_c,hA_tc_c/mcp_t_c]) - arbitrary_removal*c_t1()

    # coolant 
    dc_c1 = dT_bulkFlow(W_c,m_c_c/2,(hx_ch1_c2(t-tau_hx_c_c)+hx_ch2_c2(t-tau_hx_c_c))/2,c_c1()) + F*dT_convective([c_t1(),c_m()],c_c1(),[hA_tc_c/mcp_c_c,hA_mc_c/mcp_c_c]) - arbitrary_removal*c_c1()
    dc_c2 = dT_bulkFlow(W_c,m_c_c/2,c_c1(),c_c2()) +                                              F*dT_convective([c_t1(),c_m()],c_c2(),[hA_tc_c/mcp_c_c,hA_mc_c/mcp_c_c]) - arbitrary_removal*c_c2()

    # moderator 
    dc_m =  dT_internal(k_m,P,mcp_m_c,n()) + F*dT_convective([c_c1(),c_c2()],c_m(),[hA_mc_c/mcp_m_c,hA_mc_c/mcp_m_c])

    # FUEL-HELIUM HX1
    dhx_fh1_f1 = dT_bulkFlow(W_f,m_f_hx,c_f2(t-tau_c_hx_f),hx_fh1_f1()) + F*dT_convective([hx_fh1_t1()],hx_fh1_f1(),[hA_ft_hx/mcp_f_hx]) 
    dhx_fh1_f2 = dT_bulkFlow(W_f,m_f_hx,hx_fh1_f1(),hx_fh1_f2()) +        F*dT_convective([hx_fh1_t1()],hx_fh1_f2(),[hA_ft_hx/mcp_f_hx]) 

    # Tubes 
    dhx_fh1_t1 = dT_convective([hx_fh1_f1(),hx_fh1_f2(),hx_fh1_h1(),hx_fh1_h2()],hx_fh1_t1(),[h/mcp_t_hx for h in [hA_ft_hx,hA_ft_hx,hA_ht_hx,hA_ht_hx]])

    # Helium
    dhx_fh1_h1 = dT_bulkFlow(W_h_fh,m_h_hxfh,hx_hwf2_h2(t-tau_h_loop_f),hx_fh1_h1()) + F*dT_convective([hx_fh1_t1()],hx_fh1_h1(),[hA_ht_hx/mcp_h_hxfh])
    dhx_fh1_h2 = dT_bulkFlow(W_h_fh,m_h_hxfh,hx_fh1_h1(),hx_fh1_h2()) +         F*dT_convective([hx_fh1_t1()],hx_fh1_h2(),[hA_ht_hx/mcp_h_hxfh])

    # FUEL-HELIUM HX2 Fuel Nodes
    dhx_fh2_f1 = dT_bulkFlow(W_f,m_f_hx,c_f2(t-tau_c_hx_f), hx_fh2_f1()) + F*dT_convective([hx_fh2_t1()],hx_fh2_f1(),[hA_ft_hx/mcp_f_hx])
    dhx_fh2_f2 = dT_bulkFlow(W_f,m_f_hx,hx_fh2_f1(),hx_fh2_f2()) +         F*dT_convective([hx_fh2_t1()],hx_fh2_f2(),[hA_ft_hx/mcp_f_hx])

    # Tubes for FUEL-HELIUM HX2
    dhx_fh2_t1 = dT_convective([hx_fh2_f1(),hx_fh2_f2(),hx_fh2_h1(),hx_fh2_h2()],hx_fh2_t1(),[h/mcp_t_hx for h in [hA_ft_hx,hA_ft_hx,hA_ht_hx,hA_ht_hx]])

    # Helium for FUEL-HELIUM HX2
    dhx_fh2_h1 = dT_bulkFlow(W_h_fh,m_h_hxfh,hx_hwf1_h2(t-tau_h_btw_f),hx_fh2_h1()) + F*dT_convective([hx_fh2_t1()],hx_fh2_h1(),[hA_ht_hx/mcp_h_hxfh])
    dhx_fh2_h2 = dT_bulkFlow(W_h_fh,m_h_hxfh,hx_fh2_h1(),hx_fh2_h2()) +  F*dT_convective([hx_fh2_t1()],hx_fh2_h2(),[hA_ht_hx/mcp_h_hxfh])

    # COOLANT-HELIUM HX1
    # Fuel Nodes
    dhx_ch1_c1 = dT_bulkFlow(W_c, m_c_hx, c_c2(t-tau_c_hx_c), hx_ch1_c1()) + F*dT_convective([hx_ch1_t1()], hx_ch1_c1(), [hA_ct_hx/mcp_c_hxch])
    dhx_ch1_c2 = dT_bulkFlow(W_c, m_c_hx, hx_ch1_c1(), hx_ch1_c2()) +        F*dT_convective([hx_ch1_t1()], hx_ch1_c2(), [hA_ct_hx/mcp_c_hxch])

    # Tubes
    dhx_ch1_t1 = F*dT_convective([hx_ch1_c1(),hx_ch1_c2(),hx_ch1_h1(),hx_ch1_h2()],hx_ch1_t1(),[h/mcp_t_hxch for h in [hA_ct_hx,hA_ct_hx,hA_th_hxch,hA_th_hxch]])

    # Helium
    dhx_ch1_h1 = dT_bulkFlow(W_h_ch, m_h_hxch/2, hx_hwc1_h2(t-tau_h_loop_c), hx_ch1_h1()) + F*dT_convective([hx_ch1_t1()], hx_ch1_h1(), [hA_th_hxch/mcp_h_hxch])
    dhx_ch1_h2 = dT_bulkFlow(W_h_ch, m_h_hxch/2, hx_ch1_h1(), hx_ch1_h2()) +         F*dT_convective([hx_ch1_t1()], hx_ch1_h2(), [hA_th_hxch/mcp_h_hxch])

    # COOLANT-HELIUM HX2
    # Fuel Nodes
    dhx_ch2_c1 = dT_bulkFlow(W_c, m_c_hx, c_c2(t-tau_c_hx_c), hx_ch2_c1()) + F*dT_convective([hx_ch2_t1()], hx_ch2_c1(), [hA_ct_hx/mcp_c_hxch])
    dhx_ch2_c2 = dT_bulkFlow(W_c, m_c_hx, hx_ch2_c1(), hx_ch2_c2()) +        F*dT_convective([hx_ch2_t1()], hx_ch2_c2(), [hA_ct_hx/mcp_c_hxch])

    # Tubes
    dhx_ch2_t1 = dT_convective([hx_ch2_c1(),hx_ch2_c2(),hx_ch2_h1(),hx_ch2_h2()],hx_ch2_t1(),[h/mcp_t_hxch for h in [hA_ct_hx,hA_ct_hx,hA_th_hxch,hA_th_hxch]])

    # Helium
    dhx_ch2_h1 = dT_bulkFlow(W_h_ch, m_h_hxch/2, hx_hwc2_h2(t-tau_h_loop_c), hx_ch2_h1()) + F*dT_convective([hx_ch2_t1()], hx_ch2_h1(), [hA_th_hxch/mcp_h_hxch])
    dhx_ch2_h2 = dT_bulkFlow(W_h_ch, m_h_hxch/2, hx_ch2_h1(), hx_ch2_h2()) +         F*dT_convective([hx_ch2_t1()], hx_ch2_h2(), [hA_th_hxch/mcp_h_hxch])

    # HELIUM-WATER HX1 (FUEL LOOP)
    # Helium
    T_hhwf_h1 = dT_bulkFlow(W_h_fh, m_h_hxhw/2, hx_fh1_h2(t-tau_h_btw_f), hx_hwf1_h1()) +  F*dT_convective([hx_hwf1_t1()], hx_hwf1_h1(), [hA_ht_hxhw/mcp_h_hxhw])
    T_hhwf_h2 = dT_bulkFlow(W_h_fh, m_h_hxhw/2, hx_hwf1_h1(), hx_hwf1_h2()) + F*dT_convective([hx_hwf1_t1()], hx_hwf1_h2(), [hA_ht_hxhw/mcp_h_hxhw])

    # Tubes
    T_hhwf_t1 = F*dT_convective([hx_hwf1_h1(),hx_hwf1_h2(),hx_hwf1_w1(),hx_hwf1_w2()],hx_hwf1_t1(),[h/mcp_t_hxhw for h in [hA_ht_hxhw,hA_ht_hxhw,hA_tw_hxhw,hA_tw_hxhw]])

    # Water
    T_hhwf_w1 = dT_bulkFlow(W_hhwf_w, m_w/2, T0_hhwf_w1-hx_hwf1_w1(), hx_hwf1_w1(),dumped=True) + F*dT_convective([hx_hwf1_t1()], hx_hwf1_w1(), [hA_tw_hxhw/mcp_w])
    T_hhwf_w2 = dT_bulkFlow(W_hhwf_w, m_w/2, hx_hwf1_w1(), hx_hwf1_w2()) +                        F*dT_convective([hx_hwf1_t1()], hx_hwf1_w2(), [hA_tw_hxhw/mcp_w])

    # HELIUM-WATER HX2 (FUEL LOOP)
    # Helium
    T_hhwf2_h1 = dT_bulkFlow(W_h_fh, m_h_hxhw/2, hx_fh2_h2(t-tau_h_btw_f), hx_hwf2_h1()) +  F*dT_convective([hx_hwf2_t1()], hx_hwf2_h1(), [hA_ht_hxhw/mcp_h_hxhw])
    T_hhwf2_h2 = dT_bulkFlow(W_h_fh, m_h_hxhw/2, hx_hwf2_h1(), hx_hwf2_h2()) + F*dT_convective([hx_hwf2_t1()], hx_hwf2_h2(), [hA_ht_hxhw/mcp_h_hxhw])

    # Tubes
    T_hhwf2_t1 = F*dT_convective([hx_hwf2_h1(),hx_hwf2_h2(),hx_hwf2_w1(),hx_hwf2_w2()],hx_hwf2_t1(),[h/mcp_t_hxhw for h in [hA_ht_hxhw,hA_ht_hxhw,hA_tw_hxhw,hA_tw_hxhw]])

    # Water
    T_hhwf2_w1 = dT_bulkFlow(W_hhwf_w, m_w/2, T0_hhwf_w1-hx_hwf2_w1(), hx_hwf2_w1(),dumped=True) + F*dT_convective([hx_hwf2_t1()], hx_hwf2_w1(), [hA_tw_hxhw/mcp_w])
    T_hhwf2_w2 = dT_bulkFlow(W_hhwf_w, m_w/2, hx_hwf2_w1(), hx_hwf2_w2()) +                        F*dT_convective([hx_hwf2_t1()], hx_hwf2_w2(), [hA_tw_hxhw/mcp_w])

    # HELIUM-WATER HX1 (COOLANT LOOP)
    # Helium
    T_hhwc_h1 = dT_bulkFlow(W_h_ch, m_h_hxhwc/2, hx_ch1_h2(t-tau_h_btw_c), hx_hwc1_h1()) + F*dT_convective([hx_hwc1_t1()], hx_hwc1_h1(), [hA_ht_hxhwc/mcp_h_hxhwc])
    T_hhwc_h2 = dT_bulkFlow(W_h_ch, m_h_hxhwc/2, hx_hwc1_h1(), hx_hwc1_h2()) + F*dT_convective([hx_hwc1_t1()], hx_hwc1_h2(), [hA_ht_hxhwc/mcp_h_hxhwc])

    # Tubes
    T_hhwc_t1 = F*dT_convective([hx_hwc1_h1(),hx_hwc1_h2(),hx_hwc1_w1(),hx_hwc1_w2()],hx_hwc1_t1(),[h/mcp_t_hxhwc for h in [hA_ht_hxhwc,hA_ht_hxhwc,hA_tw_hxhwc,hA_tw_hxhwc]])

    # Water
    T_hhwc_w1 = dT_bulkFlow(W_hhwc_w, m_w_hxhwc/2, T0_hhwf_w1-hx_hwc1_w1(), hx_hwc1_w1(), dumped=True) + F*dT_convective([hx_hwc1_t1()], hx_hwc1_w1(), [hA_tw_hxhwc/mcp_w_hxhwc])
    T_hhwc_w2 = dT_bulkFlow(W_hhwc_w, m_w_hxhwc/2, hx_hwc1_w1(), hx_hwc1_w2()) + F*dT_convective([hx_hwc1_t1()], hx_hwc1_w2(), [hA_tw_hxhwc/mcp_w_hxhwc])


    # HELIUM-WATER HX2 (COOLANT LOOP)
    # Helium
    T_hhwc2_h1 = dT_bulkFlow(W_h_ch, m_h_hxhwc/2, hx_ch2_h2(t-tau_h_btw_c), hx_hwc2_h1()) + F*dT_convective([hx_hwc2_t1()], hx_hwc2_h1(), [hA_ht_hxhwc/mcp_h_hxhwc])
    T_hhwc2_h2 = dT_bulkFlow(W_h_ch, m_h_hxhwc/2, hx_hwc2_h1(), hx_hwc2_h2()) + F*dT_convective([hx_hwc2_t1()], hx_hwc2_h2(), [hA_ht_hxhwc/mcp_h_hxhwc])

    # Tubes
    T_hhwc2_t1 = F*dT_convective([hx_hwc2_h1(),hx_hwc2_h2(),hx_hwc2_w1(),hx_hwc2_w2()],hx_hwc2_t1(),[h/mcp_t_hxhwc for h in [hA_ht_hxhwc,hA_ht_hxhwc,hA_tw_hxhwc,hA_tw_hxhwc]])

    # Water
    T_hhwc2_w1 = dT_bulkFlow(W_hhwc_w, m_w_hxhwc/2, T0_hhwf_w1-hx_hwc2_w1(), hx_hwc2_w1(), dumped=True) + F*dT_convective([hx_hwc2_t1()], hx_hwc2_w1(), [hA_tw_hxhwc/mcp_w_hxhwc])
    T_hhwc2_w2 = dT_bulkFlow(W_hhwc_w, m_w_hxhwc/2, hx_hwc2_w1(), hx_hwc2_w2()) + F*dT_convective([hx_hwc2_t1()], hx_hwc2_w2(), [hA_tw_hxhwc/mcp_w_hxhwc])

    # reactivity insertion
    t0 = 300.00
    tf = 1000.00
    T = np.arange(0.0,tf,0.01)
    inserted = 4e-3
    def rho_insert(t):
        insert_duration = 0.4/0.011
        if (t<t0):
            return 0.0
        elif (t<(t0+insert_duration)):
            return ((t-t0))*(inserted/insert_duration)
        else:
            return inserted

    rho_spline = CubicHermiteSpline(n=1)
    rho_spline.from_function(rho_insert, times_of_interest = T)
    rho_ext = input(0)

    dn = ((rho()+rho_ext)-beta_t)*n()/Lam+lam[0]*C1()+lam[1]*C2()+lam[2]*C3()+lam[3]*C4()+lam[4]*C5()+lam[5]*C6()           # n (no source insertion): n()

    # dC_i/dt (precursor concentrations)
    dC1 = n()*beta_t*(beta/beta_t)[0]/Lam - lam[0]*C1() - C1()/tau_c + C1(t-tau_l)*np.exp(-lam[0]*tau_l)/tau_c                       # C1: y(27)
    dC2 = n()*beta_t*(beta/beta_t)[1]/Lam - lam[1]*C2() - C2()/tau_c + C2(t-tau_l)*np.exp(-lam[1]*tau_l)/tau_c                       # C2: y(28)
    dC3 = n()*beta_t*(beta/beta_t)[2]/Lam - lam[2]*C3() - C3()/tau_c + C3(t-tau_l)*np.exp(-lam[2]*tau_l)/tau_c                       # C3: y(29)
    dC4 = n()*beta_t*(beta/beta_t)[3]/Lam - lam[3]*C4() - C4()/tau_c + C4(t-tau_l)*np.exp(-lam[3]*tau_l)/tau_c                       # C4: y(30)
    dC5 = n()*beta_t*(beta/beta_t)[4]/Lam - lam[4]*C5() - C5()/tau_c + C5(t-tau_l)*np.exp(-lam[4]*tau_l)/tau_c                       # C5: y(31)
    dC6 = n()*beta_t*(beta/beta_t)[5]/Lam - lam[5]*C6() - C6()/tau_c + C6(t-tau_l)*np.exp(-lam[5]*tau_l)/tau_c                       # C6: y(32)

    # reactivity 
    drho = (a_f/2)*(dc_f1 + dc_f2)+(a_b)*(dc_m)+(a_c/2)*(dc_c1+dc_c2)           # rho: y(33)

    # instantiate jitcdde object
    DDE = jitcdde_input([dc_f1,dc_f2,dc_t1,dc_c1,dc_c2,dc_m,
                        dhx_fh1_f1,dhx_fh1_f2,dhx_fh1_t1,dhx_fh1_h1,dhx_fh1_h2,
                        dhx_fh2_f1,dhx_fh2_f2,dhx_fh2_t1,dhx_fh2_h1,dhx_fh2_h2,
                        dhx_ch1_c1,dhx_ch1_c2,dhx_ch1_t1,dhx_ch1_h1,dhx_ch1_h2,
                        dhx_ch2_c1,dhx_ch2_c2,dhx_ch2_t1,dhx_ch2_h1,dhx_ch2_h2,
                        T_hhwf_h1,T_hhwf_h2,T_hhwf_t1,T_hhwf_w1,T_hhwf_w2,
                        T_hhwf2_h1,T_hhwf2_h2,T_hhwf2_t1,T_hhwf2_w1,T_hhwf2_w2,
                        T_hhwc_h1,T_hhwc_h2,T_hhwc_t1,T_hhwc_w1,T_hhwc_w2,
                        T_hhwc2_h1,T_hhwc2_h2,T_hhwc2_t1,T_hhwc2_w1,T_hhwc2_w2,
                        dn,dC1,dC2,dC3,dC4,dC5,dC6,drho],
                        rho_spline)

    # set initial conditions
    DDE.constant_past([T0_c_f1,T0_c_f2,T0_c_t1,T0_c_c1,T0_c_c2,T0_c_m+50,
                    T0_hfh_f1,T0_hfh_f2,T0_hfh_t1,T0_hfh_h1,T0_hfh_h2,
                    T0_hfh_f1,T0_hfh_f2,T0_hfh_t1,T0_hfh_h1,T0_hfh_h2,
                    T0_hch_c1,T0_hch_c2,T0_hch_t1,T0_hch_h1,T0_hch_h2,
                    T0_hch_c1,T0_hch_c2,T0_hch_t1,T0_hch_h1,T0_hch_h2,
                    T0_hhwf_h1,T0_hhwf_h2,T0_hhwf_t1,T0_hhwf_w1,T0_hhwf_w2,
                    T0_hhwf_h1,T0_hhwf_h2,T0_hhwf_t1,T0_hhwf_w1,T0_hhwf_w2,
                    T0_hhwc_h1,T0_hhwc_h2,T0_hhwc_t1,T0_hhwc_w1,T0_hhwc_w2,
                    T0_hhwc_h1,T0_hhwc_h2,T0_hhwc_t1,T0_hhwc_w1,T0_hhwc_w2,
                    n_frac0,C0[0],C0[1],C0[2],C0[3],C0[4],C0[5],0.0
                        ])

    sol_jit = []
    for t_x in T:
        sol_jit.append(DDE.integrate(t_x))

    return sol_jit

def relax_ht(params):

    # Unpack parameters
    hA_ft_c, hA_tc_c, hA_mc_c, hA_ft_hx, hA_ht_hx, hA_ct_hx, hA_th_hxch, \
    hA_ht_hxhw, hA_tw_hxhw, hA_ht_hxhwc, hA_tw_hxhwc = params
    
    # CORE
    # fuel nodes
    dc_f1 = dT_bulkFlow(W_f, m_f_c/2, (hx_fh1_f2(t-tau_hx_c_f)+hx_fh2_f2(t-tau_hx_c_f))/2,c_f1()) + dT_internal(k_f1, P, mcp_f_c, n()) + F*dT_convective([c_t1()],c_f1(),[hA_ft_c/mcp_f_c]) - arbitrary_removal*c_f1()
    dc_f2 = dT_bulkFlow(W_f, m_f_c/2, c_f1(), c_f2()) +                                             dT_internal(k_f1, P, mcp_f_c, n()) + dT_convective([c_t1()],c_f2(),[hA_ft_c/mcp_f_c]) - arbitrary_removal*c_f2()

    # tubes
    dc_t1= F*dT_convective([c_f1(),c_f2(),c_c1(),c_c2()],c_t1(),[hA_ft_c/mcp_t_c,hA_ft_c/mcp_t_c,hA_tc_c/mcp_t_c,hA_tc_c/mcp_t_c]) - arbitrary_removal*c_t1()

    # coolant 
    dc_c1 = dT_bulkFlow(W_c,m_c_c/2,(hx_ch1_c2(t-tau_hx_c_c)+hx_ch2_c2(t-tau_hx_c_c))/2,c_c1()) + F*dT_convective([c_t1(),c_m()],c_c1(),[hA_tc_c/mcp_c_c,hA_mc_c/mcp_c_c]) - arbitrary_removal*c_c1()
    dc_c2 = dT_bulkFlow(W_c,m_c_c/2,c_c1(),c_c2()) +                                              F*dT_convective([c_t1(),c_m()],c_c2(),[hA_tc_c/mcp_c_c,hA_mc_c/mcp_c_c]) - arbitrary_removal*c_c2()

    # moderator 
    dc_m =  dT_internal(k_m,P,mcp_m_c,n()) + F*dT_convective([c_c1(),c_c2()],c_m(),[hA_mc_c/mcp_m_c,hA_mc_c/mcp_m_c])

    # FUEL-HELIUM HX1
    dhx_fh1_f1 = dT_bulkFlow(W_f,m_f_hx,c_f2(t-tau_c_hx_f),hx_fh1_f1()) + F*dT_convective([hx_fh1_t1()],hx_fh1_f1(),[hA_ft_hx/mcp_f_hx]) 
    dhx_fh1_f2 = dT_bulkFlow(W_f,m_f_hx,hx_fh1_f1(),hx_fh1_f2()) +        F*dT_convective([hx_fh1_t1()],hx_fh1_f2(),[hA_ft_hx/mcp_f_hx]) 

    # Tubes 
    dhx_fh1_t1 = dT_convective([hx_fh1_f1(),hx_fh1_f2(),hx_fh1_h1(),hx_fh1_h2()],hx_fh1_t1(),[h/mcp_t_hx for h in [hA_ft_hx,hA_ft_hx,hA_ht_hx,hA_ht_hx]])

    # Helium
    dhx_fh1_h1 = dT_bulkFlow(W_h_fh,m_h_hxfh,hx_hwf2_h2(t-tau_h),hx_fh1_h1()) + F*dT_convective([hx_fh1_t1()],hx_fh1_h1(),[hA_ht_hx/mcp_h_hxfh])
    dhx_fh1_h2 = dT_bulkFlow(W_h_fh,m_h_hxfh,hx_fh1_h1(),hx_fh1_h2()) +         F*dT_convective([hx_fh1_t1()],hx_fh1_h2(),[hA_ht_hx/mcp_h_hxfh])

    # FUEL-HELIUM HX2 Fuel Nodes
    dhx_fh2_f1 = dT_bulkFlow(W_f,m_f_hx,c_f2(t-tau_c_hx_f), hx_fh2_f1()) + F*dT_convective([hx_fh2_t1()],hx_fh2_f1(),[hA_ft_hx/mcp_f_hx])
    dhx_fh2_f2 = dT_bulkFlow(W_f,m_f_hx,hx_fh2_f1(),hx_fh2_f2()) +         F*dT_convective([hx_fh2_t1()],hx_fh2_f2(),[hA_ft_hx/mcp_f_hx])

    # Tubes for FUEL-HELIUM HX2
    dhx_fh2_t1 = dT_convective([hx_fh2_f1(),hx_fh2_f2(),hx_fh2_h1(),hx_fh2_h2()],hx_fh2_t1(),[h/mcp_t_hx for h in [hA_ft_hx,hA_ft_hx,hA_ht_hx,hA_ht_hx]])

    # Helium for FUEL-HELIUM HX2
    dhx_fh2_h1 = dT_bulkFlow(W_h_fh,m_h_hxfh,hx_hwf1_h2(),hx_fh2_h1()) + F*dT_convective([hx_fh2_t1()],hx_fh2_h1(),[hA_ht_hx/mcp_h_hxfh])
    dhx_fh2_h2 = dT_bulkFlow(W_h_fh,m_h_hxfh,hx_fh2_h1(),hx_fh2_h2()) +  F*dT_convective([hx_fh2_t1()],hx_fh2_h2(),[hA_ht_hx/mcp_h_hxfh])

    # COOLANT-HELIUM HX1
    # Fuel Nodes
    dhx_ch1_c1 = dT_bulkFlow(W_c, m_c_hx, c_c2(t-tau_c_hx_f), hx_ch1_c1()) + F*dT_convective([hx_ch1_t1()], hx_ch1_c1(), [hA_ct_hx/mcp_c_hxch])
    dhx_ch1_c2 = dT_bulkFlow(W_c, m_c_hx, hx_ch1_c1(), hx_ch1_c2()) +        F*dT_convective([hx_ch1_t1()], hx_ch1_c2(), [hA_ct_hx/mcp_c_hxch])

    # Tubes
    dhx_ch1_t1 = F*dT_convective([hx_ch1_c1(),hx_ch1_c2(),hx_ch1_h1(),hx_ch1_h2()],hx_ch1_t1(),[h/mcp_t_hxch for h in [hA_ct_hx,hA_ct_hx,hA_th_hxch,hA_th_hxch]])

    # Helium
    dhx_ch1_h1 = dT_bulkFlow(W_h_ch, m_h_hxch/2, hx_hwc1_h2(t-tau_h), hx_ch1_h1()) + F*dT_convective([hx_ch1_t1()], hx_ch1_h1(), [hA_th_hxch/mcp_h_hxch])
    dhx_ch1_h2 = dT_bulkFlow(W_h_ch, m_h_hxch/2, hx_ch1_h1(), hx_ch1_h2()) +         F*dT_convective([hx_ch1_t1()], hx_ch1_h2(), [hA_th_hxch/mcp_h_hxch])

    # COOLANT-HELIUM HX2
    # Fuel Nodes
    dhx_ch2_c1 = dT_bulkFlow(W_c, m_c_hx, c_c2(t-tau_c_hx_f), hx_ch2_c1()) + F*dT_convective([hx_ch2_t1()], hx_ch2_c1(), [hA_ct_hx/mcp_c_hxch])
    dhx_ch2_c2 = dT_bulkFlow(W_c, m_c_hx, hx_ch2_c1(), hx_ch2_c2()) +        F*dT_convective([hx_ch2_t1()], hx_ch2_c2(), [hA_ct_hx/mcp_c_hxch])

    # Tubes
    dhx_ch2_t1 = dT_convective([hx_ch2_c1(),hx_ch2_c2(),hx_ch2_h1(),hx_ch2_h2()],hx_ch2_t1(),[h/mcp_t_hxch for h in [hA_ct_hx,hA_ct_hx,hA_th_hxch,hA_th_hxch]])

    # Helium
    dhx_ch2_h1 = dT_bulkFlow(W_h_ch, m_h_hxch/2, hx_hwc2_h2(t-tau_h), hx_ch2_h1()) + F*dT_convective([hx_ch2_t1()], hx_ch2_h1(), [hA_th_hxch/mcp_h_hxch])
    dhx_ch2_h2 = dT_bulkFlow(W_h_ch, m_h_hxch/2, hx_ch2_h1(), hx_ch2_h2()) +         F*dT_convective([hx_ch2_t1()], hx_ch2_h2(), [hA_th_hxch/mcp_h_hxch])

    # HELIUM-WATER HX1 (FUEL LOOP)
    # Helium
    T_hhwf_h1 = dT_bulkFlow(W_h_fh, m_h_hxhw/2, hx_fh1_h2(), hx_hwf1_h1()) +  F*dT_convective([hx_hwf1_t1()], hx_hwf1_h1(), [hA_ht_hxhw/mcp_h_hxhw])
    T_hhwf_h2 = dT_bulkFlow(W_h_fh, m_h_hxhw/2, hx_hwf1_h1(), hx_hwf1_h2()) + F*dT_convective([hx_hwf1_t1()], hx_hwf1_h2(), [hA_ht_hxhw/mcp_h_hxhw])

    # Tubes
    T_hhwf_t1 = F*dT_convective([hx_hwf1_h1(),hx_hwf1_h2(),hx_hwf1_w1(),hx_hwf1_w2()],hx_hwf1_t1(),[h/mcp_t_hxhw for h in [hA_ht_hxhw,hA_ht_hxhw,hA_tw_hxhw,hA_tw_hxhw]])

    # Water
    T_hhwf_w1 = dT_bulkFlow(W_hhwf_w, m_w/2, T0_hhwf_w1-hx_hwf1_w1(), hx_hwf1_w1(),dumped=True) + F*dT_convective([hx_hwf1_t1()], hx_hwf1_w1(), [hA_tw_hxhw/mcp_w])
    T_hhwf_w2 = dT_bulkFlow(W_hhwf_w, m_w/2, hx_hwf1_w1(), hx_hwf1_w2()) +                        F*dT_convective([hx_hwf1_t1()], hx_hwf1_w2(), [hA_tw_hxhw/mcp_w])

    # HELIUM-WATER HX2 (FUEL LOOP)
    # Helium
    T_hhwf2_h1 = dT_bulkFlow(W_h_fh, m_h_hxhw/2, hx_fh2_h2(), hx_hwf2_h1()) +  F*dT_convective([hx_hwf2_t1()], hx_hwf2_h1(), [hA_ht_hxhw/mcp_h_hxhw])
    T_hhwf2_h2 = dT_bulkFlow(W_h_fh, m_h_hxhw/2, hx_hwf2_h1(), hx_hwf2_h2()) + F*dT_convective([hx_hwf2_t1()], hx_hwf2_h2(), [hA_ht_hxhw/mcp_h_hxhw])

    # Tubes
    T_hhwf2_t1 = F*dT_convective([hx_hwf2_h1(),hx_hwf2_h2(),hx_hwf2_w1(),hx_hwf2_w2()],hx_hwf2_t1(),[h/mcp_t_hxhw for h in [hA_ht_hxhw,hA_ht_hxhw,hA_tw_hxhw,hA_tw_hxhw]])

    # Water
    T_hhwf2_w1 = dT_bulkFlow(W_hhwf_w, m_w/2, T0_hhwf_w1-hx_hwf2_w1(), hx_hwf2_w1(),dumped=True) + F*dT_convective([hx_hwf2_t1()], hx_hwf2_w1(), [hA_tw_hxhw/mcp_w])
    T_hhwf2_w2 = dT_bulkFlow(W_hhwf_w, m_w/2, hx_hwf2_w1(), hx_hwf2_w2()) +                        F*dT_convective([hx_hwf2_t1()], hx_hwf2_w2(), [hA_tw_hxhw/mcp_w])

    # HELIUM-WATER HX1 (COOLANT LOOP)
    # Helium
    T_hhwc_h1 = dT_bulkFlow(W_h_ch, m_h_hxhwc/2, hx_ch1_h2(), hx_hwc1_h1()) + F*dT_convective([hx_hwc1_t1()], hx_hwc1_h1(), [hA_ht_hxhwc/mcp_h_hxhwc])
    T_hhwc_h2 = dT_bulkFlow(W_h_ch, m_h_hxhwc/2, hx_hwc1_h1(), hx_hwc1_h2()) + F*dT_convective([hx_hwc1_t1()], hx_hwc1_h2(), [hA_ht_hxhwc/mcp_h_hxhwc])

    # Tubes
    T_hhwc_t1 = F*dT_convective([hx_hwc1_h1(),hx_hwc1_h2(),hx_hwc1_w1(),hx_hwc1_w2()],hx_hwc1_t1(),[h/mcp_t_hxhwc for h in [hA_ht_hxhwc,hA_ht_hxhwc,hA_tw_hxhwc,hA_tw_hxhwc]])

    # Water
    T_hhwc_w1 = dT_bulkFlow(W_hhwc_w, m_w_hxhwc/2, T0_hhwf_w1-hx_hwc1_w1(), hx_hwc1_w1(), dumped=True) + F*dT_convective([hx_hwc1_t1()], hx_hwc1_w1(), [hA_tw_hxhwc/mcp_w_hxhwc])
    T_hhwc_w2 = dT_bulkFlow(W_hhwc_w, m_w_hxhwc/2, hx_hwc1_w1(), hx_hwc1_w2()) + F*dT_convective([hx_hwc1_t1()], hx_hwc1_w2(), [hA_tw_hxhwc/mcp_w_hxhwc])


    # HELIUM-WATER HX2 (COOLANT LOOP)
    # Helium
    T_hhwc2_h1 = dT_bulkFlow(W_h_ch, m_h_hxhwc/2, hx_ch2_h2(), hx_hwc2_h1()) + F*dT_convective([hx_hwc2_t1()], hx_hwc2_h1(), [hA_ht_hxhwc/mcp_h_hxhwc])
    T_hhwc2_h2 = dT_bulkFlow(W_h_ch, m_h_hxhwc/2, hx_hwc2_h1(), hx_hwc2_h2()) + F*dT_convective([hx_hwc2_t1()], hx_hwc2_h2(), [hA_ht_hxhwc/mcp_h_hxhwc])

    # Tubes
    T_hhwc2_t1 = F*dT_convective([hx_hwc2_h1(),hx_hwc2_h2(),hx_hwc2_w1(),hx_hwc2_w2()],hx_hwc2_t1(),[h/mcp_t_hxhwc for h in [hA_ht_hxhwc,hA_ht_hxhwc,hA_tw_hxhwc,hA_tw_hxhwc]])

    # Water
    T_hhwc2_w1 = dT_bulkFlow(W_hhwc_w, m_w_hxhwc/2, T0_hhwf_w1-hx_hwc2_w1(), hx_hwc2_w1(), dumped=True) + F*dT_convective([hx_hwc2_t1()], hx_hwc2_w1(), [hA_tw_hxhwc/mcp_w_hxhwc])
    T_hhwc2_w2 = dT_bulkFlow(W_hhwc_w, m_w_hxhwc/2, hx_hwc2_w1(), hx_hwc2_w2()) + F*dT_convective([hx_hwc2_t1()], hx_hwc2_w2(), [hA_tw_hxhwc/mcp_w_hxhwc])

    # reactivity insertion
    t0 = 300.00
    tf = 1000.00
    T = np.arange(0.0,tf,0.01)
    inserted = 4e-3
    def rho_insert(t):
        insert_duration = 0.4/0.011
        if (t<t0):
            return 0.0
        elif (t<(t0+insert_duration)):
            return ((t-t0))*(inserted/insert_duration)
        else:
            return inserted

    rho_spline = CubicHermiteSpline(n=1)
    rho_spline.from_function(rho_insert, times_of_interest = T)

    rho_ext = input(0)

    dn = ((rho()+rho_ext)-beta_t)*n()/Lam+lam[0]*C1()+lam[1]*C2()+lam[2]*C3()+lam[3]*C4()+lam[4]*C5()+lam[5]*C6()           # n (no source insertion): n()

    # dC_i/dt (precursor concentrations)
    dC1 = n()*beta_t*(beta/beta_t)[0]/Lam - lam[0]*C1() - C1()/tau_c + C1(t-tau_l)*np.exp(-lam[0]*tau_l)/tau_c                       # C1: y(27)
    dC2 = n()*beta_t*(beta/beta_t)[1]/Lam - lam[1]*C2() - C2()/tau_c + C2(t-tau_l)*np.exp(-lam[1]*tau_l)/tau_c                       # C2: y(28)
    dC3 = n()*beta_t*(beta/beta_t)[2]/Lam - lam[2]*C3() - C3()/tau_c + C3(t-tau_l)*np.exp(-lam[2]*tau_l)/tau_c                       # C3: y(29)
    dC4 = n()*beta_t*(beta/beta_t)[3]/Lam - lam[3]*C4() - C4()/tau_c + C4(t-tau_l)*np.exp(-lam[3]*tau_l)/tau_c                       # C4: y(30)
    dC5 = n()*beta_t*(beta/beta_t)[4]/Lam - lam[4]*C5() - C5()/tau_c + C5(t-tau_l)*np.exp(-lam[4]*tau_l)/tau_c                       # C5: y(31)
    dC6 = n()*beta_t*(beta/beta_t)[5]/Lam - lam[5]*C6() - C6()/tau_c + C6(t-tau_l)*np.exp(-lam[5]*tau_l)/tau_c                       # C6: y(32)

    # reactivity 
    drho = (a_f/2)*(dc_f1 + dc_f2)+(a_b)*(dc_m)+(a_c/2)*(dc_c1+dc_c2)           # rho: y(33)

    # instantiate jitcdde object
    DDE = jitcdde_input([dc_f1,dc_f2,dc_t1,dc_c1,dc_c2,dc_m,
                        dhx_fh1_f1,dhx_fh1_f2,dhx_fh1_t1,dhx_fh1_h1,dhx_fh1_h2,
                        dhx_fh2_f1,dhx_fh2_f2,dhx_fh2_t1,dhx_fh2_h1,dhx_fh2_h2,
                        dhx_ch1_c1,dhx_ch1_c2,dhx_ch1_t1,dhx_ch1_h1,dhx_ch1_h2,
                        dhx_ch2_c1,dhx_ch2_c2,dhx_ch2_t1,dhx_ch2_h1,dhx_ch2_h2,
                        T_hhwf_h1,T_hhwf_h2,T_hhwf_t1,T_hhwf_w1,T_hhwf_w2,
                        T_hhwf2_h1,T_hhwf2_h2,T_hhwf2_t1,T_hhwf2_w1,T_hhwf2_w2,
                        T_hhwc_h1,T_hhwc_h2,T_hhwc_t1,T_hhwc_w1,T_hhwc_w2,
                        T_hhwc2_h1,T_hhwc2_h2,T_hhwc2_t1,T_hhwc2_w1,T_hhwc2_w2,
                        dn,dC1,dC2,dC3,dC4,dC5,dC6,drho],
                        rho_spline)

    # set initial conditions
    DDE.constant_past([T0_c_f1,T0_c_f2,T0_c_t1,T0_c_c1,T0_c_c2,T0_c_m+50,
                    T0_hfh_f1,T0_hfh_f2,T0_hfh_t1,T0_hfh_h1,T0_hfh_h2,
                    T0_hfh_f1,T0_hfh_f2,T0_hfh_t1,T0_hfh_h1,T0_hfh_h2,
                    T0_hch_c1,T0_hch_c2,T0_hch_t1,T0_hch_h1,T0_hch_h2,
                    T0_hch_c1,T0_hch_c2,T0_hch_t1,T0_hch_h1,T0_hch_h2,
                    T0_hhwf_h1,T0_hhwf_h2,T0_hhwf_t1,T0_hhwf_w1,T0_hhwf_w2,
                    T0_hhwf_h1,T0_hhwf_h2,T0_hhwf_t1,T0_hhwf_w1,T0_hhwf_w2,
                    T0_hhwc_h1,T0_hhwc_h2,T0_hhwc_t1,T0_hhwc_w1,T0_hhwc_w2,
                    T0_hhwc_h1,T0_hhwc_h2,T0_hhwc_t1,T0_hhwc_w1,T0_hhwc_w2,
                    n_frac0,C0[0],C0[1],C0[2],C0[3],C0[4],C0[5],0.0
                        ])

    sol_jit = []
    for t_x in T:
        sol_jit.append(DDE.integrate(t_x))

    return sol_jit
    inserted = 4e-3
    def rho_insert(t):
        insert_duration = 0.4/0.011
        if (t<t0):
            return 0.0
        elif (t<(t0+insert_duration)):
            return ((t-t0))*(inserted/insert_duration)
        else:
            return inserted

    rho_spline = CubicHermiteSpline(n=1)
    rho_spline.from_function(rho_insert, times_of_interest = T)

    rho_ext = input(0)

    dn = ((rho()+rho_ext)-beta_t)*n()/Lam+lam[0]*C1()+lam[1]*C2()+lam[2]*C3()+lam[3]*C4()+lam[4]*C5()+lam[5]*C6()           # n (no source insertion): n()

    # dC_i/dt (precursor concentrations)
    dC1 = n()*beta_t*(beta/beta_t)[0]/Lam - lam[0]*C1() - C1()/tau_c + C1(t-tau_l)*np.exp(-lam[0]*tau_l)/tau_c                       # C1: y(27)
    dC2 = n()*beta_t*(beta/beta_t)[1]/Lam - lam[1]*C2() - C2()/tau_c + C2(t-tau_l)*np.exp(-lam[1]*tau_l)/tau_c                       # C2: y(28)
    dC3 = n()*beta_t*(beta/beta_t)[2]/Lam - lam[2]*C3() - C3()/tau_c + C3(t-tau_l)*np.exp(-lam[2]*tau_l)/tau_c                       # C3: y(29)
    dC4 = n()*beta_t*(beta/beta_t)[3]/Lam - lam[3]*C4() - C4()/tau_c + C4(t-tau_l)*np.exp(-lam[3]*tau_l)/tau_c                       # C4: y(30)
    dC5 = n()*beta_t*(beta/beta_t)[4]/Lam - lam[4]*C5() - C5()/tau_c + C5(t-tau_l)*np.exp(-lam[4]*tau_l)/tau_c                       # C5: y(31)
    dC6 = n()*beta_t*(beta/beta_t)[5]/Lam - lam[5]*C6() - C6()/tau_c + C6(t-tau_l)*np.exp(-lam[5]*tau_l)/tau_c                       # C6: y(32)

    # reactivity 
    drho = (a_f/2)*(dc_f1 + dc_f2)+(a_b)*(dc_m)+(a_c/2)*(dc_c1+dc_c2)           # rho: y(33)

    # instantiate jitcdde object
    DDE = jitcdde_input([dc_f1,dc_f2,dc_t1,dc_c1,dc_c2,dc_m,
                        dhx_fh1_f1,dhx_fh1_f2,dhx_fh1_t1,dhx_fh1_h1,dhx_fh1_h2,
                        dhx_fh2_f1,dhx_fh2_f2,dhx_fh2_t1,dhx_fh2_h1,dhx_fh2_h2,
                        dhx_ch1_c1,dhx_ch1_c2,dhx_ch1_t1,dhx_ch1_h1,dhx_ch1_h2,
                        dhx_ch2_c1,dhx_ch2_c2,dhx_ch2_t1,dhx_ch2_h1,dhx_ch2_h2,
                        T_hhwf_h1,T_hhwf_h2,T_hhwf_t1,T_hhwf_w1,T_hhwf_w2,
                        T_hhwf2_h1,T_hhwf2_h2,T_hhwf2_t1,T_hhwf2_w1,T_hhwf2_w2,
                        T_hhwc_h1,T_hhwc_h2,T_hhwc_t1,T_hhwc_w1,T_hhwc_w2,
                        T_hhwc2_h1,T_hhwc2_h2,T_hhwc2_t1,T_hhwc2_w1,T_hhwc2_w2,
                        dn,dC1,dC2,dC3,dC4,dC5,dC6,drho],
                        rho_spline)

    # set initial conditions
    DDE.constant_past([T0_c_f1,T0_c_f2,T0_c_t1,T0_c_c1,T0_c_c2,T0_c_m+50,
                    T0_hfh_f1,T0_hfh_f2,T0_hfh_t1,T0_hfh_h1,T0_hfh_h2,
                    T0_hfh_f1,T0_hfh_f2,T0_hfh_t1,T0_hfh_h1,T0_hfh_h2,
                    T0_hch_c1,T0_hch_c2,T0_hch_t1,T0_hch_h1,T0_hch_h2,
                    T0_hch_c1,T0_hch_c2,T0_hch_t1,T0_hch_h1,T0_hch_h2,
                    T0_hhwf_h1,T0_hhwf_h2,T0_hhwf_t1,T0_hhwf_w1,T0_hhwf_w2,
                    T0_hhwf_h1,T0_hhwf_h2,T0_hhwf_t1,T0_hhwf_w1,T0_hhwf_w2,
                    T0_hhwc_h1,T0_hhwc_h2,T0_hhwc_t1,T0_hhwc_w1,T0_hhwc_w2,
                    T0_hhwc_h1,T0_hhwc_h2,T0_hhwc_t1,T0_hhwc_w1,T0_hhwc_w2,
                    n_frac0,C0[0],C0[1],C0[2],C0[3],C0[4],C0[5],0.0
                        ])

    sol_jit = []
    for t_x in T:
        sol_jit.append(DDE.integrate(t_x))

    return sol_jit
    inserted = 4e-3
    def rho_insert(t):
        insert_duration = 0.4/0.011
        if (t<t0):
            return 0.0
        elif (t<(t0+insert_duration)):
            return ((t-t0))*(inserted/insert_duration)
        else:
            return inserted

    rho_spline = CubicHermiteSpline(n=1)
    rho_spline.from_function(rho_insert, times_of_interest = T)

    rho_ext = input(0)

    dn = ((rho()+rho_ext)-beta_t)*n()/Lam+lam[0]*C1()+lam[1]*C2()+lam[2]*C3()+lam[3]*C4()+lam[4]*C5()+lam[5]*C6()           # n (no source insertion): n()

    # dC_i/dt (precursor concentrations)
    dC1 = n()*beta_t*(beta/beta_t)[0]/Lam - lam[0]*C1() - C1()/tau_c + C1(t-tau_l)*np.exp(-lam[0]*tau_l)/tau_c                       # C1: y(27)
    dC2 = n()*beta_t*(beta/beta_t)[1]/Lam - lam[1]*C2() - C2()/tau_c + C2(t-tau_l)*np.exp(-lam[1]*tau_l)/tau_c                       # C2: y(28)
    dC3 = n()*beta_t*(beta/beta_t)[2]/Lam - lam[2]*C3() - C3()/tau_c + C3(t-tau_l)*np.exp(-lam[2]*tau_l)/tau_c                       # C3: y(29)
    dC4 = n()*beta_t*(beta/beta_t)[3]/Lam - lam[3]*C4() - C4()/tau_c + C4(t-tau_l)*np.exp(-lam[3]*tau_l)/tau_c                       # C4: y(30)
    dC5 = n()*beta_t*(beta/beta_t)[4]/Lam - lam[4]*C5() - C5()/tau_c + C5(t-tau_l)*np.exp(-lam[4]*tau_l)/tau_c                       # C5: y(31)
    dC6 = n()*beta_t*(beta/beta_t)[5]/Lam - lam[5]*C6() - C6()/tau_c + C6(t-tau_l)*np.exp(-lam[5]*tau_l)/tau_c                       # C6: y(32)

    # reactivity 
    drho = (a_f/2)*(dc_f1 + dc_f2)+(a_b)*(dc_m)+(a_c/2)*(dc_c1+dc_c2)           # rho: y(33)

    # instantiate jitcdde object
    DDE = jitcdde_input([dc_f1,dc_f2,dc_t1,dc_c1,dc_c2,dc_m,
                        dhx_fh1_f1,dhx_fh1_f2,dhx_fh1_t1,dhx_fh1_h1,dhx_fh1_h2,
                        dhx_fh2_f1,dhx_fh2_f2,dhx_fh2_t1,dhx_fh2_h1,dhx_fh2_h2,
                        dhx_ch1_c1,dhx_ch1_c2,dhx_ch1_t1,dhx_ch1_h1,dhx_ch1_h2,
                        dhx_ch2_c1,dhx_ch2_c2,dhx_ch2_t1,dhx_ch2_h1,dhx_ch2_h2,
                        T_hhwf_h1,T_hhwf_h2,T_hhwf_t1,T_hhwf_w1,T_hhwf_w2,
                        T_hhwf2_h1,T_hhwf2_h2,T_hhwf2_t1,T_hhwf2_w1,T_hhwf2_w2,
                        T_hhwc_h1,T_hhwc_h2,T_hhwc_t1,T_hhwc_w1,T_hhwc_w2,
                        T_hhwc2_h1,T_hhwc2_h2,T_hhwc2_t1,T_hhwc2_w1,T_hhwc2_w2,
                        dn,dC1,dC2,dC3,dC4,dC5,dC6,drho],
                        rho_spline)

    # set initial conditions
    DDE.constant_past([T0_c_f1,T0_c_f2,T0_c_t1,T0_c_c1,T0_c_c2,T0_c_m+50,
                    T0_hfh_f1,T0_hfh_f2,T0_hfh_t1,T0_hfh_h1,T0_hfh_h2,
                    T0_hfh_f1,T0_hfh_f2,T0_hfh_t1,T0_hfh_h1,T0_hfh_h2,
                    T0_hch_c1,T0_hch_c2,T0_hch_t1,T0_hch_h1,T0_hch_h2,
                    T0_hch_c1,T0_hch_c2,T0_hch_t1,T0_hch_h1,T0_hch_h2,
                    T0_hhwf_h1,T0_hhwf_h2,T0_hhwf_t1,T0_hhwf_w1,T0_hhwf_w2,
                    T0_hhwf_h1,T0_hhwf_h2,T0_hhwf_t1,T0_hhwf_w1,T0_hhwf_w2,
                    T0_hhwc_h1,T0_hhwc_h2,T0_hhwc_t1,T0_hhwc_w1,T0_hhwc_w2,
                    T0_hhwc_h1,T0_hhwc_h2,T0_hhwc_t1,T0_hhwc_w1,T0_hhwc_w2,
                    n_frac0,C0[0],C0[1],C0[2],C0[3],C0[4],C0[5],0.0
                        ])

    sol_jit = []
    for t_x in T:
        sol_jit.append(DDE.integrate(t_x))

    return sol_jit
    inserted = 4e-3
    def rho_insert(t):
        insert_duration = 0.4/0.011
        if (t<t0):
            return 0.0
        elif (t<(t0+insert_duration)):
            return ((t-t0))*(inserted/insert_duration)
        else:
            return inserted

    rho_spline = CubicHermiteSpline(n=1)
    rho_spline.from_function(rho_insert, times_of_interest = T)

    rho_ext = input(0)

    dn = ((rho()+rho_ext)-beta_t)*n()/Lam+lam[0]*C1()+lam[1]*C2()+lam[2]*C3()+lam[3]*C4()+lam[4]*C5()+lam[5]*C6()           # n (no source insertion): n()

    # dC_i/dt (precursor concentrations)
    dC1 = n()*beta_t*(beta/beta_t)[0]/Lam - lam[0]*C1() - C1()/tau_c + C1(t-tau_l)*np.exp(-lam[0]*tau_l)/tau_c                       # C1: y(27)
    dC2 = n()*beta_t*(beta/beta_t)[1]/Lam - lam[1]*C2() - C2()/tau_c + C2(t-tau_l)*np.exp(-lam[1]*tau_l)/tau_c                       # C2: y(28)
    dC3 = n()*beta_t*(beta/beta_t)[2]/Lam - lam[2]*C3() - C3()/tau_c + C3(t-tau_l)*np.exp(-lam[2]*tau_l)/tau_c                       # C3: y(29)
    dC4 = n()*beta_t*(beta/beta_t)[3]/Lam - lam[3]*C4() - C4()/tau_c + C4(t-tau_l)*np.exp(-lam[3]*tau_l)/tau_c                       # C4: y(30)
    dC5 = n()*beta_t*(beta/beta_t)[4]/Lam - lam[4]*C5() - C5()/tau_c + C5(t-tau_l)*np.exp(-lam[4]*tau_l)/tau_c                       # C5: y(31)
    dC6 = n()*beta_t*(beta/beta_t)[5]/Lam - lam[5]*C6() - C6()/tau_c + C6(t-tau_l)*np.exp(-lam[5]*tau_l)/tau_c                       # C6: y(32)

    # reactivity 
    drho = (a_f/2)*(dc_f1 + dc_f2)+(a_b)*(dc_m)+(a_c/2)*(dc_c1+dc_c2)           # rho: y(33)

    # instantiate jitcdde object
    DDE = jitcdde_input([dc_f1,dc_f2,dc_t1,dc_c1,dc_c2,dc_m,
                        dhx_fh1_f1,dhx_fh1_f2,dhx_fh1_t1,dhx_fh1_h1,dhx_fh1_h2,
                        dhx_fh2_f1,dhx_fh2_f2,dhx_fh2_t1,dhx_fh2_h1,dhx_fh2_h2,
                        dhx_ch1_c1,dhx_ch1_c2,dhx_ch1_t1,dhx_ch1_h1,dhx_ch1_h2,
                        dhx_ch2_c1,dhx_ch2_c2,dhx_ch2_t1,dhx_ch2_h1,dhx_ch2_h2,
                        T_hhwf_h1,T_hhwf_h2,T_hhwf_t1,T_hhwf_w1,T_hhwf_w2,
                        T_hhwf2_h1,T_hhwf2_h2,T_hhwf2_t1,T_hhwf2_w1,T_hhwf2_w2,
                        T_hhwc_h1,T_hhwc_h2,T_hhwc_t1,T_hhwc_w1,T_hhwc_w2,
                        T_hhwc2_h1,T_hhwc2_h2,T_hhwc2_t1,T_hhwc2_w1,T_hhwc2_w2,
                        dn,dC1,dC2,dC3,dC4,dC5,dC6,drho],
                        rho_spline)

    # set initial conditions
    DDE.constant_past([T0_c_f1,T0_c_f2,T0_c_t1,T0_c_c1,T0_c_c2,T0_c_m+50,
                    T0_hfh_f1,T0_hfh_f2,T0_hfh_t1,T0_hfh_h1,T0_hfh_h2,
                    T0_hfh_f1,T0_hfh_f2,T0_hfh_t1,T0_hfh_h1,T0_hfh_h2,
                    T0_hch_c1,T0_hch_c2,T0_hch_t1,T0_hch_h1,T0_hch_h2,
                    T0_hch_c1,T0_hch_c2,T0_hch_t1,T0_hch_h1,T0_hch_h2,
                    T0_hhwf_h1,T0_hhwf_h2,T0_hhwf_t1,T0_hhwf_w1,T0_hhwf_w2,
                    T0_hhwf_h1,T0_hhwf_h2,T0_hhwf_t1,T0_hhwf_w1,T0_hhwf_w2,
                    T0_hhwc_h1,T0_hhwc_h2,T0_hhwc_t1,T0_hhwc_w1,T0_hhwc_w2,
                    T0_hhwc_h1,T0_hhwc_h2,T0_hhwc_t1,T0_hhwc_w1,T0_hhwc_w2,
                    n_frac0,C0[0],C0[1],C0[2],C0[3],C0[4],C0[5],0.0
                        ])

    sol_jit = []
    for t_x in T:
        sol_jit.append(DDE.integrate(t_x))

    return sol_jit
    inserted = 4e-3
    def rho_insert(t):
        insert_duration = 0.4/0.011
        if (t<t0):
            return 0.0
        elif (t<(t0+insert_duration)):
            return ((t-t0))*(inserted/insert_duration)
        else:
            return inserted

    rho_spline = CubicHermiteSpline(n=1)
    rho_spline.from_function(rho_insert, times_of_interest = T)

    rho_ext = input(0)

    dn = ((rho()+rho_ext)-beta_t)*n()/Lam+lam[0]*C1()+lam[1]*C2()+lam[2]*C3()+lam[3]*C4()+lam[4]*C5()+lam[5]*C6()           # n (no source insertion): n()

    # dC_i/dt (precursor concentrations)
    dC1 = n()*beta_t*(beta/beta_t)[0]/Lam - lam[0]*C1() - C1()/tau_c + C1(t-tau_l)*np.exp(-lam[0]*tau_l)/tau_c                       # C1: y(27)
    dC2 = n()*beta_t*(beta/beta_t)[1]/Lam - lam[1]*C2() - C2()/tau_c + C2(t-tau_l)*np.exp(-lam[1]*tau_l)/tau_c                       # C2: y(28)
    dC3 = n()*beta_t*(beta/beta_t)[2]/Lam - lam[2]*C3() - C3()/tau_c + C3(t-tau_l)*np.exp(-lam[2]*tau_l)/tau_c                       # C3: y(29)
    dC4 = n()*beta_t*(beta/beta_t)[3]/Lam - lam[3]*C4() - C4()/tau_c + C4(t-tau_l)*np.exp(-lam[3]*tau_l)/tau_c                       # C4: y(30)
    dC5 = n()*beta_t*(beta/beta_t)[4]/Lam - lam[4]*C5() - C5()/tau_c + C5(t-tau_l)*np.exp(-lam[4]*tau_l)/tau_c                       # C5: y(31)
    dC6 = n()*beta_t*(beta/beta_t)[5]/Lam - lam[5]*C6() - C6()/tau_c + C6(t-tau_l)*np.exp(-lam[5]*tau_l)/tau_c                       # C6: y(32)

    # reactivity 
    drho = (a_f/2)*(dc_f1 + dc_f2)+(a_b)*(dc_m)+(a_c/2)*(dc_c1+dc_c2)           # rho: y(33)

    # instantiate jitcdde object
    DDE = jitcdde_input([dc_f1,dc_f2,dc_t1,dc_c1,dc_c2,dc_m,
                        dhx_fh1_f1,dhx_fh1_f2,dhx_fh1_t1,dhx_fh1_h1,dhx_fh1_h2,
                        dhx_fh2_f1,dhx_fh2_f2,dhx_fh2_t1,dhx_fh2_h1,dhx_fh2_h2,
                        dhx_ch1_c1,dhx_ch1_c2,dhx_ch1_t1,dhx_ch1_h1,dhx_ch1_h2,
                        dhx_ch2_c1,dhx_ch2_c2,dhx_ch2_t1,dhx_ch2_h1,dhx_ch2_h2,
                        T_hhwf_h1,T_hhwf_h2,T_hhwf_t1,T_hhwf_w1,T_hhwf_w2,
                        T_hhwf2_h1,T_hhwf2_h2,T_hhwf2_t1,T_hhwf2_w1,T_hhwf2_w2,
                        T_hhwc_h1,T_hhwc_h2,T_hhwc_t1,T_hhwc_w1,T_hhwc_w2,
                        T_hhwc2_h1,T_hhwc2_h2,T_hhwc2_t1,T_hhwc2_w1,T_hhwc2_w2,
                        dn,dC1,dC2,dC3,dC4,dC5,dC6,drho],
                        rho_spline)

    # set initial conditions
    DDE.constant_past([T0_c_f1,T0_c_f2,T0_c_t1,T0_c_c1,T0_c_c2,T0_c_m+50,
                    T0_hfh_f1,T0_hfh_f2,T0_hfh_t1,T0_hfh_h1,T0_hfh_h2,
                    T0_hfh_f1,T0_hfh_f2,T0_hfh_t1,T0_hfh_h1,T0_hfh_h2,
                    T0_hch_c1,T0_hch_c2,T0_hch_t1,T0_hch_h1,T0_hch_h2,
                    T0_hch_c1,T0_hch_c2,T0_hch_t1,T0_hch_h1,T0_hch_h2,
                    T0_hhwf_h1,T0_hhwf_h2,T0_hhwf_t1,T0_hhwf_w1,T0_hhwf_w2,
                    T0_hhwf_h1,T0_hhwf_h2,T0_hhwf_t1,T0_hhwf_w1,T0_hhwf_w2,
                    T0_hhwc_h1,T0_hhwc_h2,T0_hhwc_t1,T0_hhwc_w1,T0_hhwc_w2,
                    T0_hhwc_h1,T0_hhwc_h2,T0_hhwc_t1,T0_hhwc_w1,T0_hhwc_w2,
                    n_frac0,C0[0],C0[1],C0[2],C0[3],C0[4],C0[5],0.0
                        ])

    sol_jit = []
    for t_x in T:
        sol_jit.append(DDE.integrate(t_x))

    return sol_jit

def relax_all(params):


    # Unpack parameters
    a_f, a_b, a_c, inserted, beta_t, tau_hx_c_f, tau_hx_c_c, tau_c_hx_f,      \
    tau_h_loop_f, tau_h_btw_f, tau_c_hx_c,tau_h_loop_c, tau_h_btw_c, hA_ft_c, \
    hA_tc_c, hA_mc_c, hA_ft_hx, hA_ht_hx, hA_ct_hx, hA_th_hxch, hA_ht_hxhw,   \
    hA_tw_hxhw, hA_ht_hxhwc, hA_tw_hxhwc, Lam, lam0, lam1, lam2, lam3, \
    lam4, lam5 = params
    
    tau_l = tau_c_hx_f + tau_hx_c_f
    
    # CORE
    # fuel nodes
    dc_f1 = dT_bulkFlow(W_f, m_f_c/2, (hx_fh1_f2(t-tau_hx_c_f)+hx_fh2_f2(t-tau_hx_c_f))/2,c_f1()) + dT_internal(k_f1, P, mcp_f_c, n()) + F*dT_convective([c_t1()],c_f1(),[hA_ft_c/mcp_f_c]) - arbitrary_removal*c_f1()
    dc_f2 = dT_bulkFlow(W_f, m_f_c/2, c_f1(), c_f2()) +                                             dT_internal(k_f1, P, mcp_f_c, n()) + dT_convective([c_t1()],c_f2(),[hA_ft_c/mcp_f_c]) - arbitrary_removal*c_f2()

    # tubes
    dc_t1= F*dT_convective([c_f1(),c_f2(),c_c1(),c_c2()],c_t1(),[hA_ft_c/mcp_t_c,hA_ft_c/mcp_t_c,hA_tc_c/mcp_t_c,hA_tc_c/mcp_t_c]) - arbitrary_removal*c_t1()

    # coolant 
    dc_c1 = dT_bulkFlow(W_c,m_c_c/2,(hx_ch1_c2(t-tau_hx_c_c)+hx_ch2_c2(t-tau_hx_c_c))/2,c_c1()) + F*dT_convective([c_t1(),c_m1()],c_c1(),[hA_tc_c/mcp_c_c,hA_mc_c/mcp_c_c]) - arbitrary_removal*c_c1()
    dc_c2 = dT_bulkFlow(W_c,m_c_c/2,c_c1(),c_c2()) +                                              F*dT_convective([c_t1(),c_m1()],c_c2(),[hA_tc_c/mcp_c_c,hA_mc_c/mcp_c_c]) - arbitrary_removal*c_c2()

    # moderator 
    dc_m =  dT_internal(k_m,P,mcp_m_c,n()) + F*dT_convective([c_c1(),c_c2()],c_m1(),[hA_mc_c/mcp_m_c,hA_mc_c/mcp_m_c])

    # FUEL-HELIUM HX1
    dhx_fh1_f1 = dT_bulkFlow(W_f,m_f_hx,c_f2(t-tau_c_hx_f),hx_fh1_f1()) + F*dT_convective([hx_fh1_t1()],hx_fh1_f1(),[hA_ft_hx/mcp_f_hx]) 
    dhx_fh1_f2 = dT_bulkFlow(W_f,m_f_hx,hx_fh1_f1(),hx_fh1_f2()) +        F*dT_convective([hx_fh1_t1()],hx_fh1_f2(),[hA_ft_hx/mcp_f_hx]) 

    # Tubes 
    dhx_fh1_t1 = dT_convective([hx_fh1_f1(),hx_fh1_f2(),hx_fh1_h1(),hx_fh1_h2()],hx_fh1_t1(),[h/mcp_t_hx for h in [hA_ft_hx,hA_ft_hx,hA_ht_hx,hA_ht_hx]])

    # Helium
    dhx_fh1_h1 = dT_bulkFlow(W_h_fh,m_h_hxfh,hx_hwf2_h2(t-tau_h_loop_f),hx_fh1_h1()) + F*dT_convective([hx_fh1_t1()],hx_fh1_h1(),[hA_ht_hx/mcp_h_hxfh])
    dhx_fh1_h2 = dT_bulkFlow(W_h_fh,m_h_hxfh,hx_fh1_h1(),hx_fh1_h2()) +         F*dT_convective([hx_fh1_t1()],hx_fh1_h2(),[hA_ht_hx/mcp_h_hxfh])

    # FUEL-HELIUM HX2 Fuel Nodes
    dhx_fh2_f1 = dT_bulkFlow(W_f,m_f_hx,c_f2(t-tau_c_hx_f), hx_fh2_f1()) + F*dT_convective([hx_fh2_t1()],hx_fh2_f1(),[hA_ft_hx/mcp_f_hx])
    dhx_fh2_f2 = dT_bulkFlow(W_f,m_f_hx,hx_fh2_f1(),hx_fh2_f2()) +         F*dT_convective([hx_fh2_t1()],hx_fh2_f2(),[hA_ft_hx/mcp_f_hx])

    # Tubes for FUEL-HELIUM HX2
    dhx_fh2_t1 = dT_convective([hx_fh2_f1(),hx_fh2_f2(),hx_fh2_h1(),hx_fh2_h2()],hx_fh2_t1(),[h/mcp_t_hx for h in [hA_ft_hx,hA_ft_hx,hA_ht_hx,hA_ht_hx]])

    # Helium for FUEL-HELIUM HX2
    dhx_fh2_h1 = dT_bulkFlow(W_h_fh,m_h_hxfh,hx_hwf1_h2(t-tau_h_btw_f),hx_fh2_h1()) + F*dT_convective([hx_fh2_t1()],hx_fh2_h1(),[hA_ht_hx/mcp_h_hxfh])
    dhx_fh2_h2 = dT_bulkFlow(W_h_fh,m_h_hxfh,hx_fh2_h1(),hx_fh2_h2()) +  F*dT_convective([hx_fh2_t1()],hx_fh2_h2(),[hA_ht_hx/mcp_h_hxfh])

    # COOLANT-HELIUM HX1
    # Fuel Nodes
    dhx_ch1_c1 = dT_bulkFlow(W_c, m_c_hx, c_c2(t-tau_c_hx_c), hx_ch1_c1()) + F*dT_convective([hx_ch1_t1()], hx_ch1_c1(), [hA_ct_hx/mcp_c_hxch])
    dhx_ch1_c2 = dT_bulkFlow(W_c, m_c_hx, hx_ch1_c1(), hx_ch1_c2()) +        F*dT_convective([hx_ch1_t1()], hx_ch1_c2(), [hA_ct_hx/mcp_c_hxch])

    # Tubes
    dhx_ch1_t1 = F*dT_convective([hx_ch1_c1(),hx_ch1_c2(),hx_ch1_h1(),hx_ch1_h2()],hx_ch1_t1(),[h/mcp_t_hxch for h in [hA_ct_hx,hA_ct_hx,hA_th_hxch,hA_th_hxch]])

    # Helium
    dhx_ch1_h1 = dT_bulkFlow(W_h_ch, m_h_hxch/2, hx_hwc1_h2(t-tau_h_loop_c), hx_ch1_h1()) + F*dT_convective([hx_ch1_t1()], hx_ch1_h1(), [hA_th_hxch/mcp_h_hxch])
    dhx_ch1_h2 = dT_bulkFlow(W_h_ch, m_h_hxch/2, hx_ch1_h1(), hx_ch1_h2()) +         F*dT_convective([hx_ch1_t1()], hx_ch1_h2(), [hA_th_hxch/mcp_h_hxch])

    # COOLANT-HELIUM HX2
    # Fuel Nodes
    dhx_ch2_c1 = dT_bulkFlow(W_c, m_c_hx, c_c2(t-tau_c_hx_c), hx_ch2_c1()) + F*dT_convective([hx_ch2_t1()], hx_ch2_c1(), [hA_ct_hx/mcp_c_hxch])
    dhx_ch2_c2 = dT_bulkFlow(W_c, m_c_hx, hx_ch2_c1(), hx_ch2_c2()) +        F*dT_convective([hx_ch2_t1()], hx_ch2_c2(), [hA_ct_hx/mcp_c_hxch])

    # Tubes
    dhx_ch2_t1 = dT_convective([hx_ch2_c1(),hx_ch2_c2(),hx_ch2_h1(),hx_ch2_h2()],hx_ch2_t1(),[h/mcp_t_hxch for h in [hA_ct_hx,hA_ct_hx,hA_th_hxch,hA_th_hxch]])

    # Helium
    dhx_ch2_h1 = dT_bulkFlow(W_h_ch, m_h_hxch/2, hx_hwc2_h2(t-tau_h_loop_c), hx_ch2_h1()) + F*dT_convective([hx_ch2_t1()], hx_ch2_h1(), [hA_th_hxch/mcp_h_hxch])
    dhx_ch2_h2 = dT_bulkFlow(W_h_ch, m_h_hxch/2, hx_ch2_h1(), hx_ch2_h2()) +         F*dT_convective([hx_ch2_t1()], hx_ch2_h2(), [hA_th_hxch/mcp_h_hxch])

    # HELIUM-WATER HX1 (FUEL LOOP)
    # Helium
    T_hhwf_h1 = dT_bulkFlow(W_h_fh, m_h_hxhwf/2, hx_fh1_h2(t-tau_h_btw_f), hx_hwf1_h1()) +  F*dT_convective([hx_hwf1_t1()], hx_hwf1_h1(), [hA_ht_hxhw/mcp_h_hxhw])
    T_hhwf_h2 = dT_bulkFlow(W_h_fh, m_h_hxhwf/2, hx_hwf1_h1(), hx_hwf1_h2()) + F*dT_convective([hx_hwf1_t1()], hx_hwf1_h2(), [hA_ht_hxhw/mcp_h_hxhw])

    # Tubes
    T_hhwf_t1 = F*dT_convective([hx_hwf1_h1(),hx_hwf1_h2(),hx_hwf1_w1(),hx_hwf1_w2()],hx_hwf1_t1(),[h/mcp_t_hxhw for h in [hA_ht_hxhw,hA_ht_hxhw,hA_tw_hxhw,hA_tw_hxhw]])

    # Water
    T_hhwf_w1 = dT_bulkFlow(W_hhwf_w, m_w_hxhwf/2, T0_hhwf_w1-hx_hwf1_w1(), hx_hwf1_w1(),dumped=True) + F*dT_convective([hx_hwf1_t1()], hx_hwf1_w1(), [hA_tw_hxhw/mcp_w])
    T_hhwf_w2 = dT_bulkFlow(W_hhwf_w, m_w_hxhwf/2, hx_hwf1_w1(), hx_hwf1_w2()) +                        F*dT_convective([hx_hwf1_t1()], hx_hwf1_w2(), [hA_tw_hxhw/mcp_w])

    # HELIUM-WATER HX2 (FUEL LOOP)
    # Helium
    T_hhwf2_h1 = dT_bulkFlow(W_h_fh, m_h_hxhwf/2, hx_fh2_h2(t-tau_h_btw_f), hx_hwf2_h1()) +  F*dT_convective([hx_hwf2_t1()], hx_hwf2_h1(), [hA_ht_hxhw/mcp_h_hxhw])
    T_hhwf2_h2 = dT_bulkFlow(W_h_fh, m_h_hxhwf/2, hx_hwf2_h1(), hx_hwf2_h2()) + F*dT_convective([hx_hwf2_t1()], hx_hwf2_h2(), [hA_ht_hxhw/mcp_h_hxhw])

    # Tubes
    T_hhwf2_t1 = F*dT_convective([hx_hwf2_h1(),hx_hwf2_h2(),hx_hwf2_w1(),hx_hwf2_w2()],hx_hwf2_t1(),[h/mcp_t_hxhw for h in [hA_ht_hxhw,hA_ht_hxhw,hA_tw_hxhw,hA_tw_hxhw]])

    # Water
    T_hhwf2_w1 = dT_bulkFlow(W_hhwf_w, m_w_hxhwf/2, T0_hhwf_w1-hx_hwf2_w1(), hx_hwf2_w1(),dumped=True) + F*dT_convective([hx_hwf2_t1()], hx_hwf2_w1(), [hA_tw_hxhw/mcp_w])
    T_hhwf2_w2 = dT_bulkFlow(W_hhwf_w, m_w_hxhwf/2, hx_hwf2_w1(), hx_hwf2_w2()) +                        F*dT_convective([hx_hwf2_t1()], hx_hwf2_w2(), [hA_tw_hxhw/mcp_w])

    # HELIUM-WATER HX1 (COOLANT LOOP)
    # Helium
    T_hhwc_h1 = dT_bulkFlow(W_h_ch, m_h_hxhwc/2, hx_ch1_h2(t-tau_h_btw_c), hx_hwc1_h1()) + F*dT_convective([hx_hwc1_t1()], hx_hwc1_h1(), [hA_ht_hxhwc/mcp_h_hxhwc])
    T_hhwc_h2 = dT_bulkFlow(W_h_ch, m_h_hxhwc/2, hx_hwc1_h1(), hx_hwc1_h2()) + F*dT_convective([hx_hwc1_t1()], hx_hwc1_h2(), [hA_ht_hxhwc/mcp_h_hxhwc])

    # Tubes
    T_hhwc_t1 = F*dT_convective([hx_hwc1_h1(),hx_hwc1_h2(),hx_hwc1_w1(),hx_hwc1_w2()],hx_hwc1_t1(),[h/mcp_t_hxhwc for h in [hA_ht_hxhwc,hA_ht_hxhwc,hA_tw_hxhwc,hA_tw_hxhwc]])

    # Water
    T_hhwc_w1 = dT_bulkFlow(W_hhwc_w, m_w_hxhwc/2, T0_hhwf_w1-hx_hwc1_w1(), hx_hwc1_w1(), dumped=True) + F*dT_convective([hx_hwc1_t1()], hx_hwc1_w1(), [hA_tw_hxhwc/mcp_w_hxhwc])
    T_hhwc_w2 = dT_bulkFlow(W_hhwc_w, m_w_hxhwc/2, hx_hwc1_w1(), hx_hwc1_w2()) + F*dT_convective([hx_hwc1_t1()], hx_hwc1_w2(), [hA_tw_hxhwc/mcp_w_hxhwc])


    # HELIUM-WATER HX2 (COOLANT LOOP)
    # Helium
    T_hhwc2_h1 = dT_bulkFlow(W_h_ch, m_h_hxhwc/2, hx_ch2_h2(t-tau_h_btw_c), hx_hwc2_h1()) + F*dT_convective([hx_hwc2_t1()], hx_hwc2_h1(), [hA_ht_hxhwc/mcp_h_hxhwc])
    T_hhwc2_h2 = dT_bulkFlow(W_h_ch, m_h_hxhwc/2, hx_hwc2_h1(), hx_hwc2_h2()) + F*dT_convective([hx_hwc2_t1()], hx_hwc2_h2(), [hA_ht_hxhwc/mcp_h_hxhwc])

    # Tubes
    T_hhwc2_t1 = F*dT_convective([hx_hwc2_h1(),hx_hwc2_h2(),hx_hwc2_w1(),hx_hwc2_w2()],hx_hwc2_t1(),[h/mcp_t_hxhwc for h in [hA_ht_hxhwc,hA_ht_hxhwc,hA_tw_hxhwc,hA_tw_hxhwc]])

    # Water
    T_hhwc2_w1 = dT_bulkFlow(W_hhwc_w, m_w_hxhwc/2, T0_hhwf_w1-hx_hwc2_w1(), hx_hwc2_w1(), dumped=True) + F*dT_convective([hx_hwc2_t1()], hx_hwc2_w1(), [hA_tw_hxhwc/mcp_w_hxhwc])
    T_hhwc2_w2 = dT_bulkFlow(W_hhwc_w, m_w_hxhwc/2, hx_hwc2_w1(), hx_hwc2_w2()) + F*dT_convective([hx_hwc2_t1()], hx_hwc2_w2(), [hA_tw_hxhwc/mcp_w_hxhwc])

    # reactivity insertion
    # external input
    def rho_insert(t):
        if (t<t_ins):
            return 0.0
        elif (t<(t_ins+insert_duration)):
            return ((t-t_ins))*(inserted/insert_duration) # linear
        elif (t < t_wd):
            return inserted
        elif (t < t_wd+insert_duration):
            return -((t-t_wd))*(inserted/insert_duration) # linear
        else:
            return 0.0

    rho_spline = CubicHermiteSpline(n=1)
    rho_spline.from_function(rho_insert, times_of_interest = T)
    rho_ext = input(0)

    dn = ((rho()+rho_ext)-beta_t)*n()/Lam+lam0*C1()+lam1*C2()+lam2*C3()+lam3*C4()+lam4*C5()+lam5*C6()           # n (no source insertion): n()

    # dC_i/dt (precursor concentrations)
    dC1 = n()*beta_t*(beta/beta_t)[0]/Lam - lam0*C1() - C1()/tau_c + C1(t-tau_l)*np.exp(-lam[0]*tau_l)/tau_c                       # C1: y(27)
    dC2 = n()*beta_t*(beta/beta_t)[1]/Lam - lam1*C2() - C2()/tau_c + C2(t-tau_l)*np.exp(-lam[1]*tau_l)/tau_c                       # C2: y(28)
    dC3 = n()*beta_t*(beta/beta_t)[2]/Lam - lam2*C3() - C3()/tau_c + C3(t-tau_l)*np.exp(-lam[2]*tau_l)/tau_c                       # C3: y(29)
    dC4 = n()*beta_t*(beta/beta_t)[3]/Lam - lam3*C4() - C4()/tau_c + C4(t-tau_l)*np.exp(-lam[3]*tau_l)/tau_c                       # C4: y(30)
    dC5 = n()*beta_t*(beta/beta_t)[4]/Lam - lam4*C5() - C5()/tau_c + C5(t-tau_l)*np.exp(-lam[4]*tau_l)/tau_c                       # C5: y(31)
    dC6 = n()*beta_t*(beta/beta_t)[5]/Lam - lam5*C6() - C6()/tau_c + C6(t-tau_l)*np.exp(-lam[5]*tau_l)/tau_c                       # C6: y(32)

    # reactivity 
    drho = (a_f/2)*(dc_f1 + dc_f2)+(a_b)*(dc_m)+(a_c/2)*(dc_c1+dc_c2)           # rho: y(33)

    # instantiate jitcdde object
    DDE = jitcdde_input([dc_f1,dc_f2,dc_t1,dc_c1,dc_c2,dc_m,dn,dC1,dC2,dC3,dC4,dC5,dC6,drho,
                        dhx_fh1_f1,dhx_fh1_f2,dhx_fh1_t1,dhx_fh1_h1,dhx_fh1_h2,
                        dhx_fh2_f1,dhx_fh2_f2,dhx_fh2_t1,dhx_fh2_h1,dhx_fh2_h2,
                        dhx_ch1_c1,dhx_ch1_c2,dhx_ch1_t1,dhx_ch1_h1,dhx_ch1_h2,
                        dhx_ch2_c1,dhx_ch2_c2,dhx_ch2_t1,dhx_ch2_h1,dhx_ch2_h2,
                        T_hhwf_h1,T_hhwf_h2,T_hhwf_t1,T_hhwf_w1,T_hhwf_w2,
                        T_hhwf2_h1,T_hhwf2_h2,T_hhwf2_t1,T_hhwf2_w1,T_hhwf2_w2,
                        T_hhwc_h1,T_hhwc_h2,T_hhwc_t1,T_hhwc_w1,T_hhwc_w2,
                        T_hhwc2_h1,T_hhwc2_h2,T_hhwc2_t1,T_hhwc2_w1,T_hhwc2_w2],
                        rho_spline)

    # set initial conditions
    DDE.constant_past([T0_c_f1,T0_c_f2,T0_c_t1,T0_c_c1,T0_c_c2,T0_c_m+50,n_frac0,C0[0],C0[1],C0[2],C0[3],C0[4],C0[5],0.0,
                    T0_hfh_f1,T0_hfh_f2,T0_hfh_t1,T0_hfh_h1,T0_hfh_h2,
                    T0_hfh_f1,T0_hfh_f2,T0_hfh_t1,T0_hfh_h1,T0_hfh_h2,
                    T0_hch_c1,T0_hch_c2,T0_hch_t1,T0_hch_h1,T0_hch_h2,
                    T0_hch_c1,T0_hch_c2,T0_hch_t1,T0_hch_h1,T0_hch_h2,
                    T0_hhwf_h1,T0_hhwf_h2,T0_hhwf_t1,T0_hhwf_w1,T0_hhwf_w2,
                    T0_hhwf_h1,T0_hhwf_h2,T0_hhwf_t1,T0_hhwf_w1,T0_hhwf_w2,
                    T0_hhwc_h1,T0_hhwc_h2,T0_hhwc_t1,T0_hhwc_w1,T0_hhwc_w2,
                    T0_hhwc_h1,T0_hhwc_h2,T0_hhwc_t1,T0_hhwc_w1,T0_hhwc_w2
                        ])

    sol_jit = []
    for t_x in T:
        sol_jit.append(DDE.integrate(t_x))

    return sol_jit

def relax_temp(params):

    # Unpack parameters
    a_f, a_b, a_c, inserted, beta_t, hA_ft_c, hA_tc_c, hA_mc_c, hA_ft_hx, \
    hA_ht_hx, hA_ct_hx, hA_th_hxch, hA_ht_hxhw, hA_tw_hxhw, hA_ht_hxhwc, \
    hA_tw_hxhwc, scp_f, scp_c = params


    # ARE system        
    ARE = System()

    # CORE NODES
    c_f1 = Node(m = m_f_c/2, scp = scp_f, W = W_f, y0 = T0_c_f1)
    c_f2 = Node(m = m_f_c/2, scp = scp_f, W = W_f, y0 = T0_c_f2)
    c_t1 = Node(m = m_t, scp = scp_t, y0 = T0_c_t1)
    c_c1 = Node(m = m_c_c/2, scp = scp_c, W = W_c, y0 = T0_c_c1)
    c_c2 = Node(m = m_c_c/2, scp = scp_c, W = W_c, y0 = T0_c_c2) 
    c_m1 = Node(m = m_m_c, scp = scp_m, y0 = T0_c_m+50)
    n = Node(y0 = n_frac0)
    C1 = Node(y0 = C0[0])
    C2 = Node(y0 = C0[1])
    C3 = Node(y0 = C0[2])
    C4 = Node(y0 = C0[3])
    C5 = Node(y0 = C0[4])
    C6 = Node(y0 = C0[5])
    rho = Node(y0 = 0.0)

    # add reactivity input
    def rho_insert(t):
        if (t<t_ins):
            return 0.0
        elif (t<(t_ins+insert_duration)):
            return ((t-t_ins))*(inserted/insert_duration) # linear
        elif (t < t_wd):
            return inserted
        elif (t < t_wd+insert_duration):
            return inserted-((t-t_wd))*(inserted/insert_duration) # linear
        else:
            return 0.0

    rho_ext = ARE.add_input(rho_insert, T)

    # FUEL-HELIUM HX1
    hx_fh1_f1 = Node(m = m_f_hx/2, scp = scp_f, W = W_f/2, y0 = T0_hfh_f1)
    hx_fh1_f2 = Node(m = m_f_hx/2, scp = scp_f, W = W_f/2, y0 = T0_hfh_f2)
    hx_fh1_t1 = Node(m = m_t_hxfh, scp = scp_t, y0 = T0_hfh_t1)
    hx_fh1_h1 = Node(m = m_h_hxfh/2, scp = scp_h, W = W_h_fh, y0 = T0_hfh_h1)
    hx_fh1_h2 = Node(m = m_h_hxfh/2, scp = scp_h, W = W_h_fh, y0 = T0_hfh_h2)

    # FUEL-HELIUM HX2
    hx_fh2_f1 = Node(m = m_f_hx/2, scp = scp_f, W = W_f/2, y0 = T0_hfh_f1)
    hx_fh2_f2 = Node(m = m_f_hx/2, scp = scp_f, W = W_f/2, y0 = T0_hfh_f2)
    hx_fh2_t1 = Node(m = m_t_hxfh, scp = scp_t, y0 = T0_hfh_t1)
    hx_fh2_h1 = Node(m = m_h_hxfh/2, scp = scp_h, W = W_h_fh, y0 = T0_hfh_h1)
    hx_fh2_h2 = Node(m = m_h_hxfh/2, scp = scp_h, W = W_h_fh, y0 = T0_hfh_h2)

    # COOLANT-HELIUM HX1
    hx_ch1_c1 = Node(m = m_c_hx/2, scp = scp_c, W = W_c/2, y0 = T0_hch_c1)
    hx_ch1_c2 = Node(m = m_c_hx/2, scp = scp_c, W = W_c/2, y0 = T0_hch_c2)
    hx_ch1_t1 = Node(m = m_t_hxch, scp = scp_t, y0 = T0_hch_t1)
    hx_ch1_h1 = Node(m = m_h_hxch/2, scp = scp_h, W = W_h_ch, y0 = T0_hch_h1)
    hx_ch1_h2 = Node(m = m_h_hxch/2, scp = scp_h, W = W_h_ch, y0 = T0_hch_h2)

    # COOLANT-HELIUM HX2
    hx_ch2_c1 = Node(m = m_c_hx/2, scp = scp_c, W = W_c/2, y0 = T0_hch_c1)
    hx_ch2_c2 = Node(m = m_c_hx/2, scp = scp_c, W = W_c/2, y0 = T0_hch_c2)
    hx_ch2_t1 = Node(m = m_t_hxch, scp = scp_t, y0 = T0_hfh_t1)
    hx_ch2_h1 = Node(m = m_h_hxch/2, scp = scp_h, W = W_h_ch, y0 = T0_hch_h1)
    hx_ch2_h2 = Node(m = m_h_hxch/2, scp = scp_h, W = W_h_ch, y0 = T0_hch_h2)

    # HELIUM-WATER HX1 (FUEL LOOP)
    hx_hwf1_h1 = Node(m = m_h_hxhwf/2, scp = scp_h, W = W_h_fh, y0 = T0_hhwf_h1)
    hx_hwf1_h2 = Node(m = m_h_hxhwf/2, scp = scp_h, W = W_h_fh, y0 = T0_hhwf_h2)
    hx_hwf1_t1 = Node(m = m_t_hxhwf, scp = scp_t, y0 = T0_hhwf_t1)
    hx_hwf1_w1 = Node(m = m_w_hxhwf/2, scp = scp_w, W = W_hhwf_w, y0 = T0_hhwf_w1)
    hx_hwf1_w2 = Node(m = m_w_hxhwf/2, scp = scp_w, W = W_hhwf_w, y0 = T0_hhwf_w2)

    # HELIUM-WATER HX2 (FUEL LOOP)
    hx_hwf2_h1 = Node(m = m_h_hxhwf/2, scp = scp_h, W = W_h_fh, y0 = T0_hhwf_h1)
    hx_hwf2_h2 = Node(m = m_h_hxhwf/2, scp = scp_h, W = W_h_fh, y0 = T0_hhwf_h2)
    hx_hwf2_t1 = Node(m = m_t_hxhwf, scp = scp_t, y0 = T0_hhwf_t1)
    hx_hwf2_w1 = Node(m = m_w_hxhwf/2, scp = scp_w, W = W_hhwf_w, y0 = T0_hhwf_w1)
    hx_hwf2_w2 = Node(m = m_w_hxhwf/2, scp = scp_w, W = W_hhwf_w, y0 = T0_hhwf_w2)

    # HELIUM-WATER HX1 (COOLANT LOOP)
    hx_hwc1_h1 = Node(m = m_h_hxhwc/2, scp = scp_h, W = W_h_ch, y0 = T0_hhwc_h1)
    hx_hwc1_h2 = Node(m = m_h_hxhwc/2, scp = scp_h, W = W_h_ch, y0 = T0_hhwc_h2)
    hx_hwc1_t1 = Node(m = m_t_hxhwc, scp = scp_t, y0 = T0_hhwf_t1)
    hx_hwc1_w1 = Node(m = m_w_hxhwc/2, scp = scp_w, W = W_hhwc_w, y0 = T0_hhwc_w1)
    hx_hwc1_w2 = Node(m = m_w_hxhwc/2, scp = scp_w, W = W_hhwc_w, y0 = T0_hhwc_w2)

    # HELIUM-WATER HX2 (COOLANT LOOP)
    hx_hwc2_h1 = Node(m = m_h_hxhwc/2, scp = scp_h, W = W_h_ch, y0 = T0_hhwc_h1)
    hx_hwc2_h2 = Node(m = m_h_hxhwc/2, scp = scp_h, W = W_h_ch, y0 = T0_hhwc_h2)
    hx_hwc2_t1 = Node(m = m_t_hxhwc, scp = scp_t, y0 = T0_hhwf_t1)
    hx_hwc2_w1 = Node(m = m_w_hxhwc/2, scp = scp_w, W = W_hhwc_w, y0 = T0_hhwc_w1)
    hx_hwc2_w2 = Node(m = m_w_hxhwc/2, scp = scp_w, W = W_hhwc_w, y0 = T0_hhwc_w2)

    ARE.add_nodes([c_f1,c_f2,c_t1,c_c1,c_c2,c_m1,n,C1,C2,C3,C4,C5,C6,rho,
                hx_fh1_f1,hx_fh1_f2,hx_fh1_t1,hx_fh1_h1,hx_fh1_h2,
                hx_fh2_f1,hx_fh2_f2,hx_fh2_t1,hx_fh2_h1,hx_fh2_h2,
                hx_ch1_c1,hx_ch1_c2,hx_ch1_t1,hx_ch1_h1,hx_ch1_h2,
                hx_ch2_c1,hx_ch2_c2,hx_ch2_t1,hx_ch2_h1,hx_ch2_h2,
                hx_hwf1_h1,hx_hwf1_h2,hx_hwf1_t1,hx_hwf1_w1,hx_hwf1_w2,
                hx_hwf2_h1,hx_hwf2_h2,hx_hwf2_t1,hx_hwf2_w1,hx_hwf2_w2,
                hx_hwc1_h1,hx_hwc1_h2,hx_hwc1_t1,hx_hwc1_w1,hx_hwc1_w2,
                hx_hwc2_h1,hx_hwc2_h2,hx_hwc2_t1,hx_hwc2_w1,hx_hwc2_w2,
                ])
    
    # CORE
    c_f1.set_dTdt_bulkFlow(source = (hx_fh1_f2.y(t-tau_hx_c_f)+hx_fh2_f2.y(t-tau_hx_c_f))/2) 
    c_f1.set_dTdt_internal(source = n.y(), k = k_f1*P/mcp_f_c)
    c_f1.set_dTdt_convective(source = [c_t1.y()], hA = [hA_ft_c/2])

    c_f2.set_dTdt_bulkFlow(source = c_f1.y()) 
    c_f2.set_dTdt_internal(source = n.y(), k = k_f2*P/mcp_f_c)
    c_f2.set_dTdt_convective(source = [c_t1.y()], hA = [hA_ft_c/2])

    c_t1.set_dTdt_convective(source = [c_f1.y(), c_f2.y(), c_c1.y(), c_c2.y()], 
                            hA = [hA_ft_c, hA_ft_c, hA_tc_c,hA_tc_c])

    c_c1.set_dTdt_bulkFlow(source = (hx_ch1_c2.y(t-tau_c_hx_f)+hx_ch2_c2.y(t-tau_c_hx_f))/2)
    c_c1.set_dTdt_convective(source = [c_t1.y(),c_m1.y()], hA = [hA_tc_c/2,hA_mc_c/2])

    c_c2.set_dTdt_bulkFlow(source = c_c1.y())
    c_c2.set_dTdt_convective(source = [c_t1.y(),c_m1.y()], hA = [hA_tc_c/2,hA_mc_c/2])

    c_m1.set_dTdt_internal(source = n.y(), k = k_m*P/mcp_m_c)
    c_m1.set_dTdt_convective(source = [c_c1.y(), c_c2.y()], hA = [hA_mc_c]*2)

    n.set_dndt(rho.y()+rho_ext, beta_t, Lam, lam, [C1.y(), C2.y(), C3.y(), C4.y(), C5.y(), C6.y()])
    C1.set_dcdt(n.y(), beta[0], Lam, lam[0], tau_c, tau_l)
    C2.set_dcdt(n.y(), beta[1], Lam, lam[1], tau_c, tau_l)
    C3.set_dcdt(n.y(), beta[2], Lam, lam[2], tau_c, tau_l)
    C4.set_dcdt(n.y(), beta[3], Lam, lam[3], tau_c, tau_l)
    C5.set_dcdt(n.y(), beta[4], Lam, lam[4], tau_c, tau_l)
    C6.set_dcdt(n.y(), beta[5], Lam, lam[5], tau_c, tau_l)
    rho.set_drdt([c_f1.dydt(),c_f2.dydt(),c_m1.dydt(),c_c1.dydt(),c_c2.dydt()],[a_f/2,a_f/2,a_b,a_c/2,a_c/2])

    # FUEL-HELIUM HX1
    hx_fh1_f1.set_dTdt_bulkFlow(source = c_f2.y(t-tau_c_hx_f))
    hx_fh1_f1.set_dTdt_convective(source = [hx_fh1_t1.y()], hA = [hA_ft_hx/2])

    hx_fh1_f2.set_dTdt_bulkFlow(source = hx_fh1_f1.y())
    hx_fh1_f2.set_dTdt_convective(source = [hx_fh1_t1.y()], hA = [hA_ft_hx/2])

    hx_fh1_t1.set_dTdt_convective(source = [hx_fh1_f1.y(),hx_fh1_f2.y(),hx_fh1_h1.y(),hx_fh1_h2.y()],
                                hA = [hA_ft_hx, hA_ft_hx, hA_ht_hx, hA_ht_hx])

    hx_fh1_h1.set_dTdt_bulkFlow(source = hx_hwf2_h2.y(t-tau_h))
    hx_fh1_h1.set_dTdt_convective(source = [hx_fh1_t1.y()], hA = [hA_ht_hx/2])

    hx_fh1_h2.set_dTdt_bulkFlow(source = hx_fh1_h1.y())
    hx_fh1_h2.set_dTdt_convective(source = [hx_fh1_t1.y()], hA = [hA_ht_hx/2])

    # FUEL-HELIUM HX2
    hx_fh2_f1.set_dTdt_bulkFlow(source = c_f2.y(t-tau_c_hx_f))
    hx_fh2_f1.set_dTdt_convective(source = [hx_fh2_t1.y()], hA = [hA_ft_hx/2])

    hx_fh2_f2.set_dTdt_bulkFlow(source = hx_fh2_f1.y())
    hx_fh2_f2.set_dTdt_convective(source = [hx_fh2_t1.y()], hA = [hA_ft_hx/2])

    hx_fh2_t1.set_dTdt_convective(source = [hx_fh2_f1.y(),hx_fh2_f2.y(),hx_fh2_h1.y(),hx_fh2_h2.y()],
                                hA = [hA_ft_hx, hA_ft_hx, hA_ht_hx, hA_ht_hx])

    hx_fh2_h1.set_dTdt_bulkFlow(source = hx_hwf1_h2.y())
    hx_fh2_h1.set_dTdt_convective(source = [hx_fh2_t1.y()], hA = [hA_ht_hx/2])

    hx_fh2_h2.set_dTdt_bulkFlow(source = hx_fh2_h1.y())
    hx_fh2_h2.set_dTdt_convective(source = [hx_fh2_t1.y()], hA = [hA_ht_hx/2])

    # COOLANT-HELIUM HX1
    hx_ch1_c1.set_dTdt_bulkFlow(source = c_c2.y(t-tau_c_hx_f))
    hx_ch1_c1.set_dTdt_convective(source = [hx_ch1_t1.y()], hA = [hA_ct_hx/2])

    hx_ch1_c2.set_dTdt_bulkFlow(source = hx_ch1_c1.y())
    hx_ch1_c2.set_dTdt_convective(source = [hx_ch1_t1.y()], hA = [hA_ct_hx/2])

    hx_ch1_t1.set_dTdt_convective(source = [hx_ch1_c1.y(),hx_ch1_c2.y(),hx_ch1_h1.y(),hx_ch1_h2.y()],
                                hA = [hA_ct_hx, hA_ct_hx, hA_th_hxch, hA_th_hxch])

    hx_ch1_h1.set_dTdt_bulkFlow(source = hx_hwc1_h2.y(t-tau_h))
    hx_ch1_h1.set_dTdt_convective(source = [hx_ch1_t1.y()], hA = [hA_th_hxch/2])

    hx_ch1_h2.set_dTdt_bulkFlow(source = hx_ch1_h1.y())
    hx_ch1_h2.set_dTdt_convective(source = [hx_ch1_t1.y()], hA = [hA_th_hxch/2])

    # COOLANT-HELIUM HX2
    hx_ch2_c1.set_dTdt_bulkFlow(source = c_c2.y(t-tau_c_hx_f))
    hx_ch2_c1.set_dTdt_convective(source = [hx_ch2_t1.y()], hA = [hA_ct_hx/2])

    hx_ch2_c2.set_dTdt_bulkFlow(source = hx_ch2_c1.y())
    hx_ch2_c2.set_dTdt_convective(source = [hx_ch2_t1.y()], hA = [hA_ct_hx/2])

    hx_ch2_t1.set_dTdt_convective(source = [hx_ch2_c1.y(),hx_ch2_c2.y(),hx_ch2_h1.y(),hx_ch2_h2.y()],
                                hA = [hA_ct_hx, hA_ct_hx, hA_th_hxch, hA_th_hxch])

    hx_ch2_h1.set_dTdt_bulkFlow(source = hx_hwc2_h2.y(t-tau_h))
    hx_ch2_h1.set_dTdt_convective(source = [hx_ch2_t1.y()], hA = [hA_th_hxch/2])

    hx_ch2_h2.set_dTdt_bulkFlow(source = hx_ch2_h1.y())
    hx_ch2_h2.set_dTdt_convective(source = [hx_ch2_t1.y()], hA = [hA_th_hxch/2])

    # HELIUM-WATER HX1 (FUEL LOOP)
    hx_hwf1_h1.set_dTdt_bulkFlow(source = hx_fh1_h2.y())
    hx_hwf1_h1.set_dTdt_convective(source = [hx_hwf1_t1.y()], hA = [hA_ht_hxhw/2])

    hx_hwf1_h2.set_dTdt_bulkFlow(source = hx_hwf1_h1.y())
    hx_hwf1_h2.set_dTdt_convective(source = [hx_hwf1_t1.y()], hA = [hA_ht_hxhw/2])

    hx_hwf1_t1.set_dTdt_convective(source = [hx_hwf1_h1.y(),hx_hwf1_h2.y(),hx_hwf1_w1.y(),hx_hwf1_w2.y()],
                                hA = [hA_ht_hxhw,hA_ht_hxhw,hA_tw_hxhw,hA_tw_hxhw])

    hx_hwf1_w1.set_dTdt_bulkFlow(source = T0_hhwf_w1)
    hx_hwf1_w1.set_dTdt_convective(source = [hx_hwf1_t1.y()], hA = [hA_tw_hxhw/2])

    hx_hwf1_w2.set_dTdt_bulkFlow(source = hx_hwf1_w1.y())
    hx_hwf1_w2.set_dTdt_convective(source = [hx_hwf1_t1.y()], hA = [hA_tw_hxhw/2])

    # HELIUM-WATER HX2 (FUEL LOOP)
    hx_hwf2_h1.set_dTdt_bulkFlow(source = hx_fh2_h2.y())
    hx_hwf2_h1.set_dTdt_convective(source = [hx_hwf2_t1.y()], hA = [hA_ht_hxhw/2])

    hx_hwf2_h2.set_dTdt_bulkFlow(source = hx_hwf2_h1.y())
    hx_hwf2_h2.set_dTdt_convective(source = [hx_hwf2_t1.y()], hA = [hA_ht_hxhw/2])

    hx_hwf2_t1.set_dTdt_convective(source = [hx_hwf2_h1.y(),hx_hwf2_h2.y(),hx_hwf2_w1.y(),hx_hwf2_w2.y()],
                                hA = [hA_ht_hxhw,hA_ht_hxhw,hA_tw_hxhw,hA_tw_hxhw])

    hx_hwf2_w1.set_dTdt_bulkFlow(source = T0_hhwf_w1)
    hx_hwf2_w1.set_dTdt_convective(source = [hx_hwf2_t1.y()], hA = [hA_tw_hxhw/2])

    hx_hwf2_w2.set_dTdt_bulkFlow(source = hx_hwf2_w1.y())
    hx_hwf2_w2.set_dTdt_convective(source = [hx_hwf2_t1.y()], hA = [hA_tw_hxhw/2])

    # HELIUM-WATER HX1 (COOLANT LOOP)
    hx_hwc1_h1.set_dTdt_bulkFlow(source = hx_ch1_h2.y())
    hx_hwc1_h1.set_dTdt_convective(source = [hx_hwc1_t1.y()], hA = [hA_ht_hxhwc/2])

    hx_hwc1_h2.set_dTdt_bulkFlow(source = hx_hwc1_h1.y())
    hx_hwc1_h2.set_dTdt_convective(source = [hx_hwc1_t1.y()], hA = [hA_ht_hxhwc/2])

    hx_hwc1_t1.set_dTdt_convective(source = [hx_hwc1_h1.y(),hx_hwc1_h2.y(),hx_hwc1_w1.y(),hx_hwc1_w2.y()],
                                hA = [hA_ht_hxhwc,hA_ht_hxhwc,hA_tw_hxhwc,hA_tw_hxhwc])

    hx_hwc1_w1.set_dTdt_bulkFlow(source = T0_hhwc_w1)
    hx_hwc1_w1.set_dTdt_convective(source = [hx_hwc1_t1.y()], hA = [hA_tw_hxhwc/2])

    hx_hwc1_w2.set_dTdt_bulkFlow(source = hx_hwc1_w1.y())
    hx_hwc1_w2.set_dTdt_convective(source = [hx_hwc1_t1.y()], hA = [hA_tw_hxhwc/2])

    # HELIUM-WATER HX2 (COOLANT LOOP)
    hx_hwc2_h1.set_dTdt_bulkFlow(source = hx_ch2_h2.y())
    hx_hwc2_h1.set_dTdt_convective(source = [hx_hwc2_t1.y()], hA = [hA_ht_hxhwc/2])

    hx_hwc2_h2.set_dTdt_bulkFlow(source = hx_hwc2_h1.y())
    hx_hwc2_h2.set_dTdt_convective(source = [hx_hwc2_t1.y()], hA = [hA_ht_hxhwc/2])

    hx_hwc2_t1.set_dTdt_convective(source = [hx_hwc2_h1.y(),hx_hwc2_h2.y(),hx_hwc2_w1.y(),hx_hwc2_w2.y()],
                                hA = [hA_ht_hxhwc,hA_ht_hxhwc,hA_tw_hxhwc,hA_tw_hxhwc])

    hx_hwc2_w1.set_dTdt_bulkFlow(source = T0_hhwc_w1)
    hx_hwc2_w1.set_dTdt_convective(source = [hx_hwc2_t1.y()], hA = [hA_tw_hxhwc/2])

    hx_hwc2_w2.set_dTdt_bulkFlow(source = hx_hwc2_w1.y())
    hx_hwc2_w2.set_dTdt_convective(source = [hx_hwc2_t1.y()], hA = [hA_tw_hxhwc/2])

    sol_jit = ARE.solve(T)

    return sol_jit

def relax_feedback_extra_nodes(params):

    # Unpack parameters
    a_f, a_b, a_c, inserted, beta_t = params

    # ARE system        
    ARE = System()

    # CORE NODES
    c_f1 = Node(m = m_f_c/2, scp = scp_f, W = W_f, y0 = T0_c_f1)
    c_f2 = Node(m = m_f_c/2, scp = scp_f, W = W_f, y0 = T0_c_f2)
    c_t1 = Node(m = m_t, scp = scp_t, y0 = T0_c_t1)
    c_c1 = Node(m = m_c_c/2, scp = scp_c, W = W_c, y0 = T0_c_c1)
    c_c2 = Node(m = m_c_c/2, scp = scp_c, W = W_c, y0 = T0_c_c2) 
    c_m1 = Node(m = m_m_c, scp = scp_m, y0 = T0_c_m+50)
    n = Node(y0 = n_frac0)
    C1 = Node(y0 = C0[0])
    C2 = Node(y0 = C0[1])
    C3 = Node(y0 = C0[2])
    C4 = Node(y0 = C0[3])
    C5 = Node(y0 = C0[4])
    C6 = Node(y0 = C0[5])
    rho = Node(y0 = 0.0)

    # add reactivity input
    def rho_insert(t):
        if (t<t_ins):
            return 0.0
        elif (t<(t_ins+insert_duration)):
            return ((t-t_ins))*(inserted/insert_duration) # linear
        elif (t < t_wd):
            return inserted
        elif (t < t_wd+insert_duration):
            return inserted-((t-t_wd))*(inserted/insert_duration) # linear
        else:
            return 0.0

    rho_ext = ARE.add_input(rho_insert, T)

    # FUEL-HELIUM HX1
    hx_fh1_f1 = Node(m = m_f_hx/2, scp = scp_f, W = W_f/2, y0 = T0_hfh_f1)
    hx_fh1_f2 = Node(m = m_f_hx/2, scp = scp_f, W = W_f/2, y0 = T0_hfh_f2)
    hx_fh1_t1 = Node(m = m_t_hxfh, scp = scp_t, y0 = T0_hfh_t1)
    hx_fh1_h1 = Node(m = m_h_hxfh/2, scp = scp_h, W = W_h_fh, y0 = T0_hfh_h1)
    hx_fh1_h2 = Node(m = m_h_hxfh/2, scp = scp_h, W = W_h_fh, y0 = T0_hfh_h2)

    # FUEL-HELIUM HX2
    hx_fh2_f1 = Node(m = m_f_hx/2, scp = scp_f, W = W_f/2, y0 = T0_hfh_f1)
    hx_fh2_f2 = Node(m = m_f_hx/2, scp = scp_f, W = W_f/2, y0 = T0_hfh_f2)
    hx_fh2_t1 = Node(m = m_t_hxfh, scp = scp_t, y0 = T0_hfh_t1)
    hx_fh2_h1 = Node(m = m_h_hxfh/2, scp = scp_h, W = W_h_fh, y0 = T0_hfh_h1)
    hx_fh2_h2 = Node(m = m_h_hxfh/2, scp = scp_h, W = W_h_fh, y0 = T0_hfh_h2)

    # COOLANT-HELIUM HX1
    hx_ch1_c1 = Node(m = m_c_hx/2, scp = scp_c, W = W_c/2, y0 = T0_hch_c1)
    hx_ch1_c2 = Node(m = m_c_hx/2, scp = scp_c, W = W_c/2, y0 = T0_hch_c2)
    hx_ch1_t1 = Node(m = m_t_hxch, scp = scp_t, y0 = T0_hch_t1)
    hx_ch1_h1 = Node(m = m_h_hxch/2, scp = scp_h, W = W_h_ch, y0 = T0_hch_h1)
    hx_ch1_h2 = Node(m = m_h_hxch/2, scp = scp_h, W = W_h_ch, y0 = T0_hch_h2)

    # COOLANT-HELIUM HX2
    hx_ch2_c1 = Node(m = m_c_hx/2, scp = scp_c, W = W_c/2, y0 = T0_hch_c1)
    hx_ch2_c2 = Node(m = m_c_hx/2, scp = scp_c, W = W_c/2, y0 = T0_hch_c2)
    hx_ch2_t1 = Node(m = m_t_hxch, scp = scp_t, y0 = T0_hfh_t1)
    hx_ch2_h1 = Node(m = m_h_hxch/2, scp = scp_h, W = W_h_ch, y0 = T0_hch_h1)
    hx_ch2_h2 = Node(m = m_h_hxch/2, scp = scp_h, W = W_h_ch, y0 = T0_hch_h2)

    # HELIUM-WATER HX1 (FUEL LOOP)
    hx_hwf1_h1 = Node(m = m_h_hxhwf/2, scp = scp_h, W = W_h_fh, y0 = T0_hhwf_h1)
    hx_hwf1_h2 = Node(m = m_h_hxhwf/2, scp = scp_h, W = W_h_fh, y0 = T0_hhwf_h2)
    hx_hwf1_t1 = Node(m = m_t_hxhwf, scp = scp_t, y0 = T0_hhwf_t1)
    hx_hwf1_w1 = Node(m = m_w_hxhwf/2, scp = scp_w, W = W_hhwf_w, y0 = T0_hhwf_w1)
    hx_hwf1_w2 = Node(m = m_w_hxhwf/2, scp = scp_w, W = W_hhwf_w, y0 = T0_hhwf_w2)

    # HELIUM-WATER HX2 (FUEL LOOP)
    hx_hwf2_h1 = Node(m = m_h_hxhwf/2, scp = scp_h, W = W_h_fh, y0 = T0_hhwf_h1)
    hx_hwf2_h2 = Node(m = m_h_hxhwf/2, scp = scp_h, W = W_h_fh, y0 = T0_hhwf_h2)
    hx_hwf2_t1 = Node(m = m_t_hxhwf, scp = scp_t, y0 = T0_hhwf_t1)
    hx_hwf2_w1 = Node(m = m_w_hxhwf/2, scp = scp_w, W = W_hhwf_w, y0 = T0_hhwf_w1)
    hx_hwf2_w2 = Node(m = m_w_hxhwf/2, scp = scp_w, W = W_hhwf_w, y0 = T0_hhwf_w2)

    # HELIUM-WATER HX1 (COOLANT LOOP)
    hx_hwc1_h1 = Node(m = m_h_hxhwc/2, scp = scp_h, W = W_h_ch, y0 = T0_hhwc_h1)
    hx_hwc1_h2 = Node(m = m_h_hxhwc/2, scp = scp_h, W = W_h_ch, y0 = T0_hhwc_h2)
    hx_hwc1_t1 = Node(m = m_t_hxhwc, scp = scp_t, y0 = T0_hhwf_t1)
    hx_hwc1_w1 = Node(m = m_w_hxhwc/2, scp = scp_w, W = W_hhwc_w, y0 = T0_hhwc_w1)
    hx_hwc1_w2 = Node(m = m_w_hxhwc/2, scp = scp_w, W = W_hhwc_w, y0 = T0_hhwc_w2)

    # HELIUM-WATER HX2 (COOLANT LOOP)
    hx_hwc2_h1 = Node(m = m_h_hxhwc/2, scp = scp_h, W = W_h_ch, y0 = T0_hhwc_h1)
    hx_hwc2_h2 = Node(m = m_h_hxhwc/2, scp = scp_h, W = W_h_ch, y0 = T0_hhwc_h2)
    hx_hwc2_t1 = Node(m = m_t_hxhwc, scp = scp_t, y0 = T0_hhwf_t1)
    hx_hwc2_w1 = Node(m = m_w_hxhwc/2, scp = scp_w, W = W_hhwc_w, y0 = T0_hhwc_w1)
    hx_hwc2_w2 = Node(m = m_w_hxhwc/2, scp = scp_w, W = W_hhwc_w, y0 = T0_hhwc_w2)

    ARE.add_nodes([c_f1,c_f2,c_t1,c_c1,c_c2,c_m1,n,C1,C2,C3,C4,C5,C6,rho,
              hx_fh1_f1,hx_fh1_f2,hx_fh1_t1,hx_fh1_h1,hx_fh1_h2,
              hx_fh2_f1,hx_fh2_f2,hx_fh2_t1,hx_fh2_h1,hx_fh2_h2,
              hx_ch1_c1,hx_ch1_c2,hx_ch1_t1,hx_ch1_h1,hx_ch1_h2,
              hx_ch2_c1,hx_ch2_c2,hx_ch2_t1,hx_ch2_h1,hx_ch2_h2,
              hx_hwf1_h1,hx_hwf1_h2,hx_hwf1_t1,hx_hwf1_w1,hx_hwf1_w2,
              hx_hwf2_h1,hx_hwf2_h2,hx_hwf2_t1,hx_hwf2_w1,hx_hwf2_w2,
              hx_hwc1_h1,hx_hwc1_h2,hx_hwc1_t1,hx_hwc1_w1,hx_hwc1_w2,
              hx_hwc2_h1,hx_hwc2_h2,hx_hwc2_t1,hx_hwc2_w1,hx_hwc2_w2,
              ])
    
    # CORE
    c_f1.set_dTdt_bulkFlow(source = (hx_fh1_f2.y(t-tau_hx_c_f)+hx_fh2_f2.y(t-tau_hx_c_f))/2) 
    c_f1.set_dTdt_internal(source = n.y(), k = k_f1*P/mcp_f_c)
    c_f1.set_dTdt_convective(source = [c_t1.y()], hA = [hA_ft_c/2])

    c_f2.set_dTdt_bulkFlow(source = c_f1.y()) 
    c_f2.set_dTdt_internal(source = n.y(), k = k_f2*P/mcp_f_c)
    c_f2.set_dTdt_convective(source = [c_t1.y()], hA = [hA_ft_c/2])

    c_t1.set_dTdt_convective(source = [c_f1.y(), c_f2.y(), c_c1.y(), c_c2.y()], 
                            hA = [hA_ft_c/2,hA_ft_c/2,hA_tc_c/2,hA_tc_c/2])

    c_c1.set_dTdt_bulkFlow(source = (hx_ch1_c2.y(t-tau_c_hx_f)+hx_ch2_c2.y(t-tau_c_hx_f))/2)
    c_c1.set_dTdt_convective(source = [c_t1.y(),c_m1.y()], hA = [hA_tc_c/2,hA_mc_c/2])

    c_c2.set_dTdt_bulkFlow(source = c_c1.y())
    c_c2.set_dTdt_convective(source = [c_t1.y(),c_m1.y()], hA = [hA_tc_c/2,hA_mc_c/2])

    c_m1.set_dTdt_internal(source = n.y(), k = k_m*P/mcp_m_c)
    c_m1.set_dTdt_convective(source = [c_c1.y(), c_c2.y()], hA = [hA_mc_c/2]*2)

    n.set_dndt(rho.y()+rho_ext, beta_t, Lam, lam, [C1.y(), C2.y(), C3.y(), C4.y(), C5.y(), C6.y()])
    C1.set_dcdt(n.y(), beta[0], Lam, lam[0], tau_c, tau_l)
    C2.set_dcdt(n.y(), beta[1], Lam, lam[1], tau_c, tau_l)
    C3.set_dcdt(n.y(), beta[2], Lam, lam[2], tau_c, tau_l)
    C4.set_dcdt(n.y(), beta[3], Lam, lam[3], tau_c, tau_l)
    C5.set_dcdt(n.y(), beta[4], Lam, lam[4], tau_c, tau_l)
    C6.set_dcdt(n.y(), beta[5], Lam, lam[5], tau_c, tau_l)
    rho.set_drdt([c_f1.dydt(),c_f2.dydt(),c_m1.dydt(),c_c1.dydt(),c_c2.dydt()],[a_f/2,a_f/2,a_b,a_c/2,a_c/2])

    # FUEL-HELIUM HX1
    hx_fh1_f1.set_dTdt_bulkFlow(source = c_f2.y(t-tau_c_hx_f))
    hx_fh1_f1.set_dTdt_convective(source = [hx_fh1_t1.y()], hA = [hA_ft_hx/2])

    hx_fh1_f2.set_dTdt_bulkFlow(source = hx_fh1_f1.y())
    hx_fh1_f2.set_dTdt_convective(source = [hx_fh1_t1.y()], hA = [hA_ft_hx/2])

    hx_fh1_t1.set_dTdt_convective(source = [hx_fh1_f1.y(),hx_fh1_f2.y(),hx_fh1_h1.y(),hx_fh1_h2.y()],
                                hA = [hA_ft_hx/2,hA_ft_hx/2,hA_ht_hx/2,hA_ht_hx/2])

    hx_fh1_h1.set_dTdt_bulkFlow(source = hx_hwf2_h2.y(t-tau_h))
    hx_fh1_h1.set_dTdt_convective(source = [hx_fh1_t1.y()], hA = [hA_ht_hx/2])

    hx_fh1_h2.set_dTdt_bulkFlow(source = hx_fh1_h1.y())
    hx_fh1_h2.set_dTdt_convective(source = [hx_fh1_t1.y()], hA = [hA_ht_hx/2])

    # FUEL-HELIUM HX2
    hx_fh2_f1.set_dTdt_bulkFlow(source = c_f2.y(t-tau_c_hx_f))
    hx_fh2_f1.set_dTdt_convective(source = [hx_fh2_t1.y()], hA = [hA_ft_hx/2])

    hx_fh2_f2.set_dTdt_bulkFlow(source = hx_fh2_f1.y())
    hx_fh2_f2.set_dTdt_convective(source = [hx_fh2_t1.y()], hA = [hA_ft_hx/2])

    hx_fh2_t1.set_dTdt_convective(source = [hx_fh2_f1.y(),hx_fh2_f2.y(),hx_fh2_h1.y(),hx_fh2_h2.y()],
                                hA = [hA_ft_hx/2,hA_ft_hx/2,hA_ht_hx/2,hA_ht_hx/2])

    hx_fh2_h1.set_dTdt_bulkFlow(source = hx_hwf1_h2.y())
    hx_fh2_h1.set_dTdt_convective(source = [hx_fh2_t1.y()], hA = [hA_ht_hx/2])

    hx_fh2_h2.set_dTdt_bulkFlow(source = hx_fh2_h1.y())
    hx_fh2_h2.set_dTdt_convective(source = [hx_fh2_t1.y()], hA = [hA_ht_hx/2])

    # COOLANT-HELIUM HX1
    hx_ch1_c1.set_dTdt_bulkFlow(source = c_c2.y(t-tau_c_hx_f))
    hx_ch1_c1.set_dTdt_convective(source = [hx_ch1_t1.y()], hA = [hA_ct_hx/2])

    hx_ch1_c2.set_dTdt_bulkFlow(source = hx_ch1_c1.y())
    hx_ch1_c2.set_dTdt_convective(source = [hx_ch1_t1.y()], hA = [hA_ct_hx/2])

    hx_ch1_t1.set_dTdt_convective(source = [hx_ch1_c1.y(),hx_ch1_c2.y(),hx_ch1_h1.y(),hx_ch1_h2.y()],
                                hA = [hA_ct_hx/2,hA_ct_hx/2,hA_th_hxch/2,hA_th_hxch/2])

    hx_ch1_h1.set_dTdt_bulkFlow(source = hx_hwc1_h2.y(t-tau_h))
    hx_ch1_h1.set_dTdt_convective(source = [hx_ch1_t1.y()], hA = [hA_th_hxch/2])

    hx_ch1_h2.set_dTdt_bulkFlow(source = hx_ch1_h1.y())
    hx_ch1_h2.set_dTdt_convective(source = [hx_ch1_t1.y()], hA = [hA_th_hxch/2])

    # COOLANT-HELIUM HX2
    hx_ch2_c1.set_dTdt_bulkFlow(source = c_c2.y(t-tau_c_hx_f))
    hx_ch2_c1.set_dTdt_convective(source = [hx_ch2_t1.y()], hA = [hA_ct_hx/2])

    hx_ch2_c2.set_dTdt_bulkFlow(source = hx_ch2_c1.y())
    hx_ch2_c2.set_dTdt_convective(source = [hx_ch2_t1.y()], hA = [hA_ct_hx/2])

    hx_ch2_t1.set_dTdt_convective(source = [hx_ch2_c1.y(),hx_ch2_c2.y(),hx_ch2_h1.y(),hx_ch2_h2.y()],
                                hA = [hA_ct_hx/2,hA_ct_hx/2,hA_th_hxch/2,hA_th_hxch/2])

    hx_ch2_h1.set_dTdt_bulkFlow(source = hx_hwc2_h2.y(t-tau_h))
    hx_ch2_h1.set_dTdt_convective(source = [hx_ch2_t1.y()], hA = [hA_th_hxch/2])

    hx_ch2_h2.set_dTdt_bulkFlow(source = hx_ch2_h1.y())
    hx_ch2_h2.set_dTdt_convective(source = [hx_ch2_t1.y()], hA = [hA_th_hxch/2])

    # HELIUM-WATER HX1 (FUEL LOOP)
    hx_hwf1_h1.set_dTdt_bulkFlow(source = hx_fh1_h2.y())
    hx_hwf1_h1.set_dTdt_convective(source = [hx_hwf1_t1.y()], hA = [hA_ht_hxhw/2])

    hx_hwf1_h2.set_dTdt_bulkFlow(source = hx_hwf1_h1.y())
    hx_hwf1_h2.set_dTdt_convective(source = [hx_hwf1_t1.y()], hA = [hA_ht_hxhw/2])

    hx_hwf1_t1.set_dTdt_convective(source = [hx_hwf1_h1.y(),hx_hwf1_h2.y(),hx_hwf1_w1.y(),hx_hwf1_w2.y()],
                                hA = [hA_ht_hxhw/2,hA_ht_hxhw/2,hA_tw_hxhw/2,hA_tw_hxhw/2])

    hx_hwf1_w1.set_dTdt_bulkFlow(source = T0_hhwf_w1)
    hx_hwf1_w1.set_dTdt_convective(source = [hx_hwf1_t1.y()], hA = [hA_tw_hxhw/2])

    hx_hwf1_w2.set_dTdt_bulkFlow(source = hx_hwf1_w1.y())
    hx_hwf1_w2.set_dTdt_convective(source = [hx_hwf1_t1.y()], hA = [hA_tw_hxhw/2])

    # HELIUM-WATER HX2 (FUEL LOOP)
    hx_hwf2_h1.set_dTdt_bulkFlow(source = hx_fh2_h2.y())
    hx_hwf2_h1.set_dTdt_convective(source = [hx_hwf2_t1.y()], hA = [hA_ht_hxhw/2])

    hx_hwf2_h2.set_dTdt_bulkFlow(source = hx_hwf2_h1.y())
    hx_hwf2_h2.set_dTdt_convective(source = [hx_hwf2_t1.y()], hA = [hA_ht_hxhw/2])

    hx_hwf2_t1.set_dTdt_convective(source = [hx_hwf2_h1.y(),hx_hwf2_h2.y(),hx_hwf2_w1.y(),hx_hwf2_w2.y()],
                                hA = [hA_ht_hxhw/2,hA_ht_hxhw/2,hA_tw_hxhw/2,hA_tw_hxhw/2])

    hx_hwf2_w1.set_dTdt_bulkFlow(source = T0_hhwf_w1)
    hx_hwf2_w1.set_dTdt_convective(source = [hx_hwf2_t1.y()], hA = [hA_tw_hxhw/2])

    hx_hwf2_w2.set_dTdt_bulkFlow(source = hx_hwf2_w1.y())
    hx_hwf2_w2.set_dTdt_convective(source = [hx_hwf2_t1.y()], hA = [hA_tw_hxhw/2])

    # HELIUM-WATER HX1 (COOLANT LOOP)
    hx_hwc1_h1.set_dTdt_bulkFlow(source = hx_ch1_h2.y())
    hx_hwc1_h1.set_dTdt_convective(source = [hx_hwc1_t1.y()], hA = [hA_ht_hxhwc/2])

    hx_hwc1_h2.set_dTdt_bulkFlow(source = hx_hwc1_h1.y())
    hx_hwc1_h2.set_dTdt_convective(source = [hx_hwc1_t1.y()], hA = [hA_ht_hxhwc/2])

    hx_hwc1_t1.set_dTdt_convective(source = [hx_hwc1_h1.y(),hx_hwc1_h2.y(),hx_hwc1_w1.y(),hx_hwc1_w2.y()],
                                hA = [hA_ht_hxhwc/2,hA_ht_hxhwc/2,hA_tw_hxhwc/2,hA_tw_hxhwc/2])

    hx_hwc1_w1.set_dTdt_bulkFlow(source = T0_hhwc_w1)
    hx_hwc1_w1.set_dTdt_convective(source = [hx_hwc1_t1.y()], hA = [hA_tw_hxhwc/2])

    hx_hwc1_w2.set_dTdt_bulkFlow(source = hx_hwc1_w1.y())
    hx_hwc1_w2.set_dTdt_convective(source = [hx_hwc1_t1.y()], hA = [hA_tw_hxhwc/2])

    # HELIUM-WATER HX2 (COOLANT LOOP)
    hx_hwc2_h1.set_dTdt_bulkFlow(source = hx_ch2_h2.y())
    hx_hwc2_h1.set_dTdt_convective(source = [hx_hwc2_t1.y()], hA = [hA_ht_hxhwc/2])

    hx_hwc2_h2.set_dTdt_bulkFlow(source = hx_hwc2_h1.y())
    hx_hwc2_h2.set_dTdt_convective(source = [hx_hwc2_t1.y()], hA = [hA_ht_hxhwc/2])

    hx_hwc2_t1.set_dTdt_convective(source = [hx_hwc2_h1.y(),hx_hwc2_h2.y(),hx_hwc2_w1.y(),hx_hwc2_w2.y()],
                                hA = [hA_ht_hxhwc/2,hA_ht_hxhwc/2,hA_tw_hxhwc/2,hA_tw_hxhwc/2])

    hx_hwc2_w1.set_dTdt_bulkFlow(source = T0_hhwc_w1)
    hx_hwc2_w1.set_dTdt_convective(source = [hx_hwc2_t1.y()], hA = [hA_tw_hxhwc/2])

    hx_hwc2_w2.set_dTdt_bulkFlow(source = hx_hwc2_w1.y())
    hx_hwc2_w2.set_dTdt_convective(source = [hx_hwc2_t1.y()], hA = [hA_tw_hxhwc/2])

    # redfine tube nodes 
    c_t1.m = m_t/2
    c_t2 = Node(m = m_t/2, scp = scp_t, y0 = T0_c_t1)

    # add second tube node to system 
    ARE.add_nodes([c_t2])

    # update dynamics
    c_f1.set_dTdt_convective(source = [c_t1.y()], hA = [hA_ft_c/2])
    c_f2.set_dTdt_convective(source = [c_t2.y()], hA = [hA_ft_c/2])
    c_t1.set_dTdt_convective(source = [c_f1.y(),c_c2.y()],hA = [hA_ft_c/2,hA_tc_c/2])
    c_t2.set_dTdt_convective(source = [c_f2.y(),c_c1.y()],hA = [hA_ft_c/2,hA_tc_c/2])
    c_c1.set_dTdt_convective(source = [c_t2.y(),c_m1.y()],hA = [hA_tc_c/2,hA_mc_c/2])
    c_c2.set_dTdt_convective(source = [c_t1.y(),c_m1.y()],hA = [hA_tc_c/2,hA_mc_c/2])

    # needs to be redeclared for changes to propogate to feedback
    rho.set_drdt([c_f1.dydt(),c_f2.dydt(),c_m1.dydt(),c_c1.dydt(),c_c2.dydt()],[a_f/2,a_f/2,a_b,a_c/2,a_c/2])

    sol_jit = ARE.solve(T)
    return sol_jit