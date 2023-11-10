% instantiate parameters 
Parameters_MSRE_U235_1region;

lags = [tau_c tau_c_hx tau_hx_c tau_hx_r tau_l tau_r_hx];
y0 = [T0_rp; T0_rs; T0_p1; T0_p2; T0_p3; T0_p4; T0_t1; T0_t2; T0_s1; T0_s2; 
          T0_s3; T0_s4; n_frac0; C0(1); C0(2); C0(3); C0(4); C0(5); C0(6); 
          T0_g1; T0_f1; T0_f2];

tspan = [0 500];
sol = dde23(@ddefun, lags, y0, tspan);

% unpack python data 
% filename = 'sim_out_500.0_8.txt'; % specify the file name
% solM = readmatrix(filename);


%yint = deval(sol,tint);

%filename = 'simout_k_500.xlsx';
%writematrix(sol.x,sol.y(13,:),filename);