%----------change parameters to colaborate with AL6082T6-------------
clear;clc;close all;
format long;
W0=1e9;
b=5;
a=0.4;                %sensitivity of sequence effect(control alp>0)
E=69.4e9;              %Young's modulus
k=500e6;                 %hardening parameter
nu=0.33;                     %poisson's ratio
y=2.98e8;         %macroscopic yield stress
stepnumber=100;
delta_alp=0.01;
lamratio=0.5; 
cycles=1;          %numerical cycles to get the mean value
cycles90=5;       %90 out of phase cycles to mean
save('AL6082T6.mat','W0','b','a','E','k','nu','y','stepnumber','cycles','delta_alp','cycles90');
load('gaussian.mat');
%--------find the best iteration lambda to get alp_m(best beta implicit)
%--------1st attempt------
delta_lam=0.1;
lamplus_range=[0.9:delta_lam:1.3]'
least_norm=zeros(size(lamplus_range));
j=1;
for  j=1:length(lamplus_range)
lamplus=lamplus_range(j);
lamminus=lamratio.*lamplus;
run('AL6082T6_bt1D_new_alp.m');%--to get first alp(p) using b(0)---
run('AL6082T6_bt1D_converge_alp.m');%--iterate alp and b to converge b
least_norm(j)=resnorm
end
[M1,I1]=min(least_norm) %find least square of lam value and index
lamplus=lamplus_range(I1) %1st best fit lambda least square
% 
% %--------2nd attempt(narrow down to next precision)------
% lamplus_range=[lamplus-delta_lam:1/3*delta_lam:lamplus+delta_lam]'
% least_norm=zeros(size(lamplus_range));
% j=1;
% for  j=1:length(lamplus_range)
% lamplus=lamplus_range(j);
% lamminus=lamratio.*lamplus;
% run('AL6082T6_bt1D_new_alp.m');
% run('AL6082T6_bt1D_converge_alp.m');
% least_norm(j)=resnorm
% end
% [M2,I2]=min(least_norm) %find least square of lam value and index
% lamplus=lamplus_range(I2);  %2nd best fit lambda least square

lamminus=lamratio.*lamplus;
W0 
b 
lamplus
lamminus
save('AL6082T6.mat','lamplus','lamminus','W0','b','alp_ben','alp_tor','a','E','k','nu','y','stepnumber','cycles','delta_alp','cycles90'); %final fitting
%---------plot------------
run('AL6082T6_bt1D_plot.m'); 

% run('AL6082T6_bt2D_err1.m')
% run('AL6082T6_bt2D_90_err1.m')
