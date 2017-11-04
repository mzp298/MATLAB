%----------change parameters to colaborate with AL6082T6-------------
clear;clc;close all;
format long;
W0=1e9;
b=1.2;
lamratio=1;
lamplus=0;
lamminus=0;
fb=1.1;
a=0.4;                %sensitivity of sequence effect(control alp>0)
E=69.4e9;              %Young's modulus
k=5e8;       %k=5e8 originally   %hardening parameter
nu=0.33;                     %poisson's ratio
y=2.98e8;         %macroscopic yield stress
stepnumber=100;
delta_alp=1e-4;
cycles=2;          %numerical cycles to get the mean value
cycles90=3;       %90 out of phase cycles to mean

lamplus_num=0.9; %manual numerical fit
save('AL6082T6.mat','W0','b','fb','delta_alp','a','lamplus_num','lamratio','lamplus','lamminus','E','k','nu','y','stepnumber','cycles','cycles90');
%------------torsion fit-------------------------
run('AL6082T6_bt1D_new_alp.m');%--to get alp Smax hydro
run('AL6082T6_bt1D_converge_alp_num.m');%--to get NF_num
run('AL6082T6_bt1D_plot.m');

%%------numerical calculation
run('AL6082T6_bt2D_err1_num.m')
run('AL6082T6_bt2D_90_err1_num.m')

run('AL6082T6_plot_allcases.m')


