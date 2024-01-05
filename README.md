# Source Codes of Real-time Derivative Pricing and Hedging with Consistent Metamodels
This repository contains the source codes used in the following paper:  
- Guangxin Jiang, L. Jeff Hong, and Shen, Haihui (2023). Real-time Derivative Pricing and Hedging with Consistent Metamodels. Submitted to _INFORMS Journal on Computing_.
## Disclaimer
**LICENSE**: Redistribution and use in source and binary forms, with or without modification, are permitted, under the terms of the BSD license.  
**WARNING**: These codes are written only for the purpose of demonstration and verification. While the correctness has been carefully checked, the quality such as standardability, clarity, generality, and efficiency has not been we considered.  
## 1 Introduction
The \src folder contains the entire MATLAB implementations of all the numerical experiments in Jiang et al. (2023). The codes are split into the following three folders.  
- \Example1 contains codes for numerical experiments in Section 6.1.  
- \Example2 contains codes for numerical experiments in Section 6.2.  
- \Example3 contains codes for numerical experiments in Section 6.3.  
## 2 Environment
The codes were written and run in MATLAB R2022a, on Windows 10 Professional 64-bit OS, with lntel(R) Xeon(R) Gold 6248R CPU\@3.00GHz. 256GB RAM
## 3 Details on implementation
In Section 6.1, Figures 1-4 are obtained by implementing the following file in \Example1\ Margin\_compare:  
-- KrigEuro_cokriging_deri.m  
-- KrigEuro_cokriging_deri_vega.m   
-- KrigEuro_cokriging_deri_theta.m  
-- KrigEuro_cokriging_deri_rho.m  
Figure 5 is obtained by implementing the following file in \Example1\IMSE  
-- IMSE_all.m  
Figure 6 is obtained by implementing the following file in \Example1\hedging\_effect  
-- hedge_var.m  
Figure 7 is obtained by implementing the following file in \Example1\hedging\_effect  
-- hedge_cost.m

In Section 6.2, Figure 8 is obtained by implementing the following file in \Example2\IMSE  
-- IMSE_VG.m  
Figure 9 is obtained by implementing the following file in \Example2\hedging\_effect  
-- hedge_var.m  

In Section 6.3, Figure 10 and Table 2 are obtained by implementing the following file in \Example3\Example3-1\  
-- run_BS.m  
Figures 11-12 and Table 3 are obtained by implementing the following file in \Example3\Example3-2\  
-- run40.m

In Section EC.9 (E-Companion), Figures EC.1 and EC.2 are obtained by implementing the following file in \EC\IMSE\_bias\_var\  
-- IMSE_all.m 

In Section EC.10 (E-Companion), Figures EC.3 and EC.7 are obtained by implementing the following file in \EC\LHS\_IMSE\  
-- IMSE_all.m 
