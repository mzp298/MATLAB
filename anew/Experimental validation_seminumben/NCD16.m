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
delta_alp=0.01;
cycles=1;          %numerical cycles to get the mean value
cycles90=5;
lamplus_num=0; %manual numerical fit
save('NCD16.mat','W0','b','fb','a','lamplus_num','lamratio','lamplus','lamminus','E','k','nu','y','stepnumber','cycles','cycles90','-append'); 
load('gaussian.mat');
%------------torsion fit-------------------------
run('NCD16_bt1D_new_alp.m');%--to get alp, Smax, hydro
run('NCD16_bt1D_converge_alp_num.m');%--identify W0, b
run('NCD16_bt1D_plot.m');

lamplus=1.4
lamminus=lamplus-1
save('NCD16.mat','lamplus','lamminus','-append'); 

run('NCD16_b1D_m_err1.m'); 
run('NCD16_bt2D_m_err1.m')
run('NCD16_bt2D_m_90_err1.m')

% %--------find the best iteration lambda to get alp_m(best beta implicit)
% %--------1st attempt------
% delta_lam=0.1;
% lamplus_range=[0:delta_lam:0]'
% least_norm=zeros(size(lamplus_range));
% j=1;
% for  j=1:length(lamplus_range)
% lamplus=lamplus_range(j);
% lamminus=lamratio.*lamplus;
% save('NCD16.mat','lamplus','lamminus','-append');
% run('NCD16_bt1D_new_alp.m');%--to get first alp(p) using b(0)---
% run('NCD16_bt1D_converge_alp.m');%--iterate alp and b to converge b
% least_norm(j)=resnorm
% end
% [M1,I1]=min(least_norm) %find least square of lam value and index
% lamplus=lamplus_range(I1) %1st best fit lambda least square

% %--------2nd attempt(narrow down to next precision)------
% lamplus_range=[lamplus-delta_lam:1/3*delta_lam:lamplus+delta_lam]'
% least_norm=zeros(size(lamplus_range));
% j=1;
% for  j=1:length(lamplus_range)
% lamplus=lamplus_range(j);
% lamminus=lamratio.*lamplus;
% save('NCD16.mat','lamplus','lamminus','-append');
% run('NCD16_bt1D_new_alp.m');
% run('NCD16_bt1D_converge_alp.m');
% least_norm(j)=resnorm
% end
% [M2,I2]=min(least_norm) %find least square of lam value and index
% lamplus=lamplus_range(I2);  %2nd best fit lambda least square

% lamminus=lamratio.*lamplus;
% W0 
% b 
% lamplus
% lamminus
% save('NCD16.mat','lamplus','lamminus','W0','b','-append');

% clear;clc;close all;
% run('NCD16_bt1D_new_alp.m');%--to get alp Smax hydro



