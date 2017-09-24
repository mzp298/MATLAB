%----------change parameters to colaborate with SM45C STEEL(HCF part)-------------
clear;clc;close all;
format long;
W0=1e7;
b=2;
lamratio=1; 
lamplus=0.5;
lamminus=lamratio.*lamplus;
fb=1.1;
a=0.01;                %sensitivity of sequence effect(control alp>0)
E=213e9;              %Young's modulus
k=1e9;                 %hardening parameter
nu=0.3;                     %poisson's ratio
y=638e6;         %macroscopic yield stress
stepnumber=100;
delta_alp=0.01;
cycles=1;          %numerical cycles to get the mean value
cycles90=5;       %90 out of phase cycles to mean
lamplus_num=0; %manual numerical fit
save('SM45C.mat','W0','b','fb','a','lamplus_num','lamratio','lamplus','lamminus','E','k','nu','y','stepnumber','cycles','cycles90','-append'); 
%------------torsion fit-------------------------

run('SM45C_bt1D_new_alp.m');%--to get alp, Smax, hydro
run('SM45C_bt1D_converge_alp_num.m');%--identify W0, b
run('SM45C_bt1D_plot.m'); 

% run('SM45C_b1D_m_err1.m') %-identify lam+ with b1D_m
% run('SM45C_bt2D_90_err1.m')
% run('SM45C_bt2D_m_90_err1.m')

%--------find the best iteration lambda to get alp_m(best beta implicit)
% %--------1st attempt------
% delta_lam=0.2;
% lamplus_range=[1.2:delta_lam:1.2]'
% least_norm=zeros(size(lamplus_range));
% j=1;
% for  j=1:length(lamplus_range)
% lamplus=lamplus_range(j);
% lamminus=lamratio.*lamplus;
% save('SM45C.mat','lamplus','lamminus','-append');
% run('SM45C_bt1D_new_alp.m');%--to get first alp(p) using b(0)---
% run('SM45C_bt1D_converge_alp.m');%--iterate alp and b to converge b
% least_norm(j)=resnorm
% end
% [M1,I1]=min(least_norm) %find least square of lam value and index
% lamplus=lamplus_range(I1) %1st best fit lambda least square
% 
% %--------2nd attempt------
% lamplus_range=[lamplus-delta_lam:1/3*delta_lam:lamplus+delta_lam]'
% least_norm=zeros(size(lamplus_range));
% j=1;
% for  j=1:length(lamplus_range)
% lamplus=lamplus_range(j);
% lamminus=lamratio.*lamplus;
% save('SM45C.mat','lamplus','lamminus','-append');
% run('SM45C_bt1D_new_alp.m');
% run('SM45C_bt1D_converge_alp.m');
% least_norm(j)=resnorm
% end
% [M2,I2]=min(least_norm) %find least square of lam value and index
% lamplus=lamplus_range(I2);  %2nd best fit lambda least square
% lamminus=lamratio.*lamplus;
% 
% W0 
% b
% lamplus
% lamminus=lamratio.*lamplus
% save('SM45C.mat','lamplus','lamminus','W0','b','-append');

% run('SM45C_b1D_m_err1.m') %-identify lam+,b,W0 with b1D_m
% clear;clc;close all;
% run('SM45C_bt1D_new_alp.m');%--to get alp Smax hydro
% run('SM45C_bt1D_converge_alp.m');%--to get NF_num
% run('SM45C_bt1D_plot.m'); 
% run('SM45C_bt2D_90_err1.m')
% run('SM45C_bt2D_m_90_err1.m')
