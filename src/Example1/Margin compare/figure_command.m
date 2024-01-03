%% record figure command
%% Delta S0 change [80,120], K=105, r=0.02,sigma=0.2, T=1
KrigEuro_cokriging_deri(100,0.02,0.2,1,105,10,1,10000,0,0,0)
%% Vega S0=100, K=105, r=0.02,sigma change [0,0.3], T=1
KrigEuro_cokriging_deri_vega(100,0.02,10,1,105,10,1,10000,0,0,0)
%% Theta S0=100, K=105, r=0.02,sigma=0.2, T change [0,1]
KrigEuro_cokriging_deri_theta(100,0.02,0.2,1,105,10,1,10000,0,0,0)
%% Theta S0=100, K=105, r=[0.001,0.1],sigma change=0.2, T=1
KrigEuro_cokriging_deri_rho(100,0.02,0.2,1,105,10,1,10000,0,0,0)