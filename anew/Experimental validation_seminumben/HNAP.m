%----------change parameters to colaborate with 10HNAP STEEL-------------
clear;clc;close all;
format long;
W0=1e9;
b=5.3;
lamplus=0;
lamratio=0.1;
lamminus=lamratio*lamplus;
fb=1.1;
a=0.4;              %sensitivity of sequence effect(control alp>0)
E=215e9;              %Young's modulus
k=8e8;                 %hardening parameter
nu=0.29;                     %poisson's ratio
y=4.18e8;         %macroscopic yield stress
stepnumber=100;
cycles=2;          %numerical cycles to get the mean value
cycles90=5;       %90 out of phase cycles to mean
lamplus_num=1.7; %manual numerical fit
save('HNAP.mat','W0','b','lamplus_num','lamplus','lamminus','fb','a','lamratio','E','k','nu','y','stepnumber','cycles','cycles90','-append');
%------------torsion fit-------------------------
run('HNAP_bt1D_new_alp.m');
run('HNAP_bt1D_converge_alp_num.m');%--to get NF_num
run('HNAP_bt1D_plot.m'); 

run('HNAP_b1D_m_err1_num.m')

% clear;clc;close all;
% run('HNAP_bt1D_new_alp.m');
% run('HNAP_bt1D_converge_alp.m');%--to get NF_num
% run('HNAP_bt1D_plot.m'); 
%% 
% 
% %-----------bt fit-----------------------------
% delta_lam=0.1;
% lamplus_range=[0.1:delta_lam:0.1]'
% least_norm=zeros(size(lamplus_range));
% j=1;
% for  j=1:length(lamplus_range)
% lamplus=lamplus_range(j);
% lamminus=lamratio.*lamplus;
% save('HNAP.mat','lamplus','lamminus','-append');
% run('HNAP_bt1D_new_alp.m');%--to get first alp(p) using b(0)---
% run('HNAP_bt1D_converge_alp.m');%--iterate alp and b to converge b
% least_norm(j)=resnorm
% end
% [M1,I1]=min(least_norm); %find least square of lam value and index
% lamplus=lamplus_range(I1) %1st best fit lambda least square
% % %------------------------
% % % %--------2nd attempt(narrow down to next precision)------
% % lamplus_range=[lamplus-delta_lam:1/3*delta_lam:lamplus+delta_lam]'
% % least_norm=zeros(size(lamplus_range));
% % j=1;
% % for  j=1:length(lamplus_range)
% % lamplus=lamplus_range(j);
% % lamminus=lamratio.*lamplus;
% % run('HNAP_bt1D_new_alp.m');%--to get first alp(p) using b(0)---
% % run('HNAP_bt1D_converge_alp.m');%--iterate alp and b to converge b
% % least_norm(j)=resnorm
% % end
% % [M2,I2]=min(least_norm) %200step 5cycles =1.553407226713651e+02
% % lamplus=lamplus_range(I2);  %2nd best fit lambda least square
% %%------------------------
% W0 
% b 
% lamplus
% lamminus
% 
% save('HNAP.mat','lamplus','lamminus','-append'); %final fitting
save('F:\Git\MATLAB\anew\Experimental validation_lam+-\HNAP_random\HNAP_random.mat','lamplus','lamminus','W0','b','fb','a','E','k','nu','y','stepnumber');



% run('HNAP_b1D_m_err1.m')