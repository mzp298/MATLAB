clear;clc;
close all
format long e
run('Parameters_ep_a_b')
load('gaussian.mat');
%vpa(x,289);
%vpa(weight,289);
% fid = fopen('F:\Git\Cetim\ep_b_05\Acqui_CV.txt');
fid = fopen('/home/ma/MATLAB/Cetim/ep_b_05/Acqui_CV.txt');
area=2.97*9.99*1e-6; %meter square
NF=11947368;
run('Referece_scalar_iteration');
xlwrite('b_fitting.xls',nF,1,'D05');
xlwrite('b_fitting.xls',(NF-nF)/NF,1,'E05');
% p70=length(find(stress11(1:repetition)>70e6))/repetition;
% xlwrite('/home/ma/a_b1.1_a0.1_lam0.1_ari/greater_stress_proportion.xlsx',p70,1,'B11');
% p90=length(find(stress11(1:repetition)>90e6))/repetition;
% xlwrite('/home/ma/a_b1.1_a0.1_lam0.1_ari/greater_stress_proportion.xlsx',p90,1,'C11');
% p110=length(find(stress11(1:repetition)>110e6))/repetition;
% xlwrite('/home/ma/a_b1.1_a0.1_lam0.1_ari/greater_stress_proportion.xlsx',p110,1,'D11');
% p130=length(find(stress11(1:repetition)>130e6))/repetition;
% xlwrite('/home/ma/a_b1.1_a0.1_lam0.1_ari/greater_stress_proportion.xlsx',p130,1,'E11');
% p150=length(find(stress11(1:repetition)>150e6))/repetition;
% xlwrite('/home/ma/a_b1.1_a0.1_lam0.1_ari/greater_stress_proportion.xlsx',p150,1,'F11');
% p170=length(find(stress11(1:repetition)>170e6))/repetition;
% xlwrite('/home/ma/a_b1.1_a0.1_lam0.1_ari/greater_stress_proportion.xlsx',p170,1,'G11');
% p190=length(find(stress11(1:repetition)>190e6))/repetition;
% xlwrite('/home/ma/a_b1.1_a0.1_lam0.1_ari/greater_stress_proportion.xlsx',p190,1,'H11');
