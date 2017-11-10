%----------change parameters to colaborate with 10HNAP STEEL-------------
clear;clc;close all;
format long;
W0=1e9;
b=5.3;
lamplus=0;
lamratio=0.1;
lamminus=lamratio*lamplus;
fb=1.1;
a=0.01;              %sensitivity of sequence effect(control alp>0)
lamplus_num=1.7; %manual numerical fit
E=215e9;              %Young's modulus
k=8e8;                 %hardening parameter
nu=0.29;                     %poisson's ratio
y=4.18e8;         %macroscopic yield stress
stepnumber=100;
delta_alp=1e-4;
cycles=2;          %numerical cycles to get the mean value
cycles90=5;       %90 out of phase cycles to mean
save('HNAP.mat','W0','b','delta_alp','lamplus_num','lamplus','lamminus','fb','a','lamratio','E','k','nu','y','stepnumber','cycles','cycles90','-append');
%------------torsion fit-------------------------
run('HNAP_bt1D_new_alp.m');
run('HNAP_bt1D_converge_alp_num.m');%--to get NF_num
run('HNAP_bt1D_plot.m');
%% 

run('HNAP_b1D_m_err1_num.m')
run('HNAP_plot_allcases.m')

% save('F:\Git\MATLAB\anew\Experimental validation_numall\HNAP_random\HNAP_random.mat','lamplus','lamminus','W0','b','fb','a','E','k','nu','y','stepnumber');
