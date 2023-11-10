function dydt = ddefun(t,y,Z)

    Parameters_MSRE_U235_1region;

    % unpack delay variables
    T_hc4_delay = Z(12,4);
    T_cf2_delay = Z(22,2);
    T_out_rc_delay = Z(1,6);
    C1_delay = Z(14,5);
    C2_delay = Z(15,5);
    C3_delay = Z(16,5);
    C4_delay = Z(17,5);
    C5_delay = Z(18,5);
    C6_delay = Z(19,5);
    T_hf4_delay = Z(6,3);
    
    % reactivity
    rho = (a_f/2)*((-T0_f1+y(21))+(-T0_f2+y(22))) + a_g*(-T0_g1+y(20)) + rho_ext;
    
    % derivatives 
    dydt = [ 
    (W_rp/mn_rp)*(T_hc4_delay-y(1))+(hA_rpn/mcp_rpn)*(y(2)-y(1));                                    % T_out_rc: y(1)
    -((W_rs/mn_rs)+(hA_rsn/mcp_rsn))*y(2)+(hA_rsn/mcp_rsn)*y(1)+ (W_rs/mn_rs)*Trs_in;                % T_out_air: y(2)
    -((W_p/mn_p)+(hA_pn/mcp_pn))*y(3)+(hA_pn/mcp_pn)*y(7)+(W_p/mn_p)*T_cf2_delay;                    % T_hf1: y(3)
    (W_p/mn_p)*(y(3)-y(4)) + (hA_pn/mcp_pn)*(y(7)-y(3));                                             % T_hf2: y(4)
    -((W_p/mn_p)+(hA_pn/mcp_pn))*y(5)+(hA_pn/mcp_pn)*y(8)+(W_p/mn_p)*y(4);                           % T_hf3: y(5)
    (W_p/mn_p)*(y(5)-y(6)) + (hA_pn/mcp_pn)*(y(8)-y(5));                                             % T_hf4: y(6)
    (2*hA_pn/mcp_tn)*(y(3)-y(7))+(2*hA_sn/mcp_tn)*(y(11)-y(7));                                      % T_ht1: y(7)
    (2*hA_pn/mcp_tn)*(y(5)-y(8))+(2*hA_sn/mcp_tn)*(y(9)-y(8));                                       % T_ht2: y(8)
    -((W_s/mn_s)+(hA_sn/mcp_sn))*y(9)+(hA_sn/mcp_sn)*y(8)+(W_s/mn_s)*T_out_rc_delay;                 % T_hc1: y(9)
    (W_s/mn_s)*(y(9)-y(10))+(hA_sn/mcp_sn)*(y(8)-y(9));                                              % T_hc2: y(10)
    -((W_s/mn_s)+(hA_sn/mcp_sn))*y(11)+(hA_sn/mcp_sn)*y(7)+(W_s/mn_s)*y(10);                         % T_hc3: y(11)
    (W_s/mn_s)*(y(11)-y(12)) + (hA_sn/mcp_sn)*(y(7)-y(11));                                          % T_hc4: y(12)
    (rho-beta_t)*y(13)/Lam+lam(1)*y(14)+lam(2)*y(15)+lam(3)*y(16)+lam(4)*y(17)+lam(5)*y(18)+lam(6)*y(19);% n: y(13)
    y(13)*beta_vec(1)/Lam-lam(1)*y(14)-y(14)/tau_c+C1_delay*exp(-lam(1)*tau_l)/tau_c;                 % C1: y(14) 
    y(13)*beta_vec(2)/Lam-lam(2)*y(15)-y(15)/tau_c+C2_delay*exp(-lam(2)*tau_l)/tau_c;                 % C2: y(15)
    y(13)*beta_vec(3)/Lam-lam(3)*y(16)-y(16)/tau_c+C3_delay*exp(-lam(3)*tau_l)/tau_c;                 % C3: y(16)
    y(13)*beta_vec(4)/Lam-lam(4)*y(17)-y(17)/tau_c+C4_delay*exp(-lam(4)*tau_l)/tau_c;                 % C4: y(17)
    y(13)*beta_vec(5)/Lam-lam(5)*y(18)-y(18)/tau_c+C5_delay*exp(-lam(5)*tau_l)/tau_c;                 % C5: y(18)
    y(13)*beta_vec(6)/Lam-lam(6)*y(19)-y(19)/tau_c+C6_delay*exp(-lam(6)*tau_l)/tau_c;                 % C6: y(19)
    (hA_fg/mcp_g1)*(y(21)-y(20))+k_g*P*y(13)/mcp_g1;                                                 % T_cg: y(20)
    W_f/mn_f*(T_hf4_delay-y(21)) + (k_f1*P*y(13)/mcp_f1) + (hA_fg*k_1*(y(20)-y(21))/mcp_f1);         % T_cf1: y(21)   
    W_f/mn_f*(y(21)-y(22))+(k_f2*P*y(13)/mcp_f2)+(hA_fg*k_2*(y(20)-y(21))/mcp_f2)                    % T_cf2: y(22)   
    ];

end