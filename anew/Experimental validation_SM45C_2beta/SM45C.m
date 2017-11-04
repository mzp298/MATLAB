%----------change parameters to colaborate with SM45C STEEL(HCF part)-------------
clear;clc;close all;
format long;
W0=1e7;
b=2;
fb=1.1;
a=0.1;                %sensitivity of sequence effect(control alp>0)
E=213e9;              %Young's modulus
k=1e9;                 %hardening parameter
nu=0.3;                     %poisson's ratio
y=638e6;         %macroscopic yield stress
stepnumber=100;
delta_alp=1e-4;
cycles=2;          %numerical cycles to get the mean value
cycles90=5;       %90 out of phase cycles to mean
lamplus=1.4; %manual numerical fit
lamminus=0;
ff=1e6*442;
tt=1e6*311;
save('SM45C.mat','tt','ff','W0','b','delta_alp','fb','a',...
    'lamplus','lamminus','E','k','nu','y','stepnumber','cycles','cycles90'); 

%------------LCF torsion fit-------------------------
run('SM45C_bt1D_new_alp_LCF.m');%--to get alp, Smax, hydro
resnorm_limit=4e12;
ub_beta=6;
run('SM45C_bt1D_converge_alp_num.m');%--identify W0, b
save('SM45C.mat','W0','b','-append'); %--1dm in LCF regime
run('SM45C_bt1D_plot_LCF.m'); 
NFben_LCF=NFben;
NFtor_LCF=NFtor;
NFben_num_LCF=NFben_num;
NFtor_num_LCF=NFtor_num;
save('SM45C.mat','NFben_LCF','NFtor_LCF','NFben_num_LCF','NFtor_num_LCF','-append');
%%
%------------HCF torsion fit-------------------------
run('SM45C_bt1D_new_alp_HCF.m');%--to get alp, Smax, hydro
resnorm_limit=4e11;
ub_beta=50;
run('SM45C_bt1D_converge_alp_num.m');%--identify W0, b
run('SM45C_bt1D_plot_HCF.m'); 
NFben_HCF=NFben;
NFtor_HCF=NFtor;
NFben_num_HCF=NFben_num;
NFtor_num_HCF=NFtor_num;
save('SM45C.mat','NFben_HCF','NFtor_HCF','NFben_num_HCF','NFtor_num_HCF','-append');

%%
% %%------numerical calculation
run('SM45C_b1D_m_err1_num.m')
run('SM45C_bt2D_90_err1_num.m')
run('SM45C_bt2D_m_90_err1_num.m')
run('SM45C_plot_allcases.m')
