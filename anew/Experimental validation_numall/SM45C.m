%----------change parameters to colaborate with SM45C STEEL(HCF part)-------------
clear;clc;close all;
format long;
W0=1e7;
b=2;
lamratio=1; 
lamplus=0;
lamminus=lamratio.*lamplus;
fb=1.1;
a=0.4;                %sensitivity of sequence effect(control alp>0)
E=213e9;              %Young's modulus
k=1e9;                 %hardening parameter
nu=0.3;                     %poisson's ratio
y=638e6;         %macroscopic yield stress
stepnumber=100;
delta_alp=1e-4;
cycles=1;          %numerical cycles to get the mean value
cycles90=5;       %90 out of phase cycles to mean
lamplus_num=1; %manual numerical fit
save('SM45C.mat','W0','b','delta_alp','fb','a','lamplus_num','lamratio','lamplus','lamminus','E','k','nu','y','stepnumber','cycles','cycles90'); 
%------------torsion fit-------------------------
run('SM45C_bt1D_new_alp.m');%--to get alp, Smax, hydro
run('SM45C_bt1D_converge_alp_num.m');%--identify W0, b
run('SM45C_bt1D_plot.m'); 

%%
% %%------analytical formula-----
% run('SM45C_b1D_m_err1.m')
% run('SM45C_bt2D_90_err1.m')
% run('SM45C_bt2D_m_90_err1.m')

% %%------numerical calculation
run('SM45C_b1D_m_err1_num.m')
run('SM45C_bt2D_90_err1_num.m')
run('SM45C_bt2D_m_90_err1_num.m')
