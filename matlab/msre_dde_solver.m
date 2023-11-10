% instantiate parameters 
Parameters_MSRE_U235_1region;

lags = [tau_c tau_c_hx tau_hx_c tau_hx_r tau_l tau_r_hx];

tspan = [0 500];
sol = dde23(@ddefun, lags, @history, tspan);


yint = deval(sol,tint);

filename = 'simout_k_500.xlsx';
writematrix(sol.x,sol.y(13,:),filename);