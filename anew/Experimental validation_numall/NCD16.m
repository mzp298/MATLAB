%----------change parameters to colaborate with NCD16 STEEL-------------
clear;clc;close all;
format long;
W0=3e8;
b=1.2;
lamratio=1; 
lamplus=0.5;
lamminus=lamratio.*lamplus;
fb=1.1;
a=0.4;                %sensitivity of sequence effect(control alp>0)
E=191e9;              %Young's modulus
k=1e9;                 %hardening parameter
nu=0.29;                     %poisson's ratio
y=1.08e9;         %macroscopic yield stress
stepnumber=100;
delta_alp=1e-4;
cycles=2;          %numerical cycles to get the mean value
cycles90=5;
lamplus_num=0.55; %manual numerical fit
save('NCD16.mat','W0','b','delta_alp','fb','a','lamplus_num','lamratio','lamplus','lamminus','E','k','nu','y','stepnumber','cycles','cycles90','-append'); 
load('gaussian.mat');
%------------torsion fit-------------------------
run('NCD16_bt1D_new_alp.m');%--to get alp, Smax, hydro
run('NCD16_bt1D_converge_alp_num.m');%--identify W0, b
run('NCD16_bt1D_plot.m');
%% 

% %%------analytical formula-----
% run('NCD16_b1D_m_err1.m'); 
% run('NCD16_bt2D_m_err1.m')
% run('NCD16_bt2D_m_90_err1.m')

%%------numerical calculation
run('NCD16_b1D_m_err1_num.m'); 
run('NCD16_bt2D_m_err1_num.m')
run('NCD16_bt2D_m_90_err1_num.m')