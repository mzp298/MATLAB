clear;clc;
close all
format long e
run('Parameters_ep_a_b')
load('gaussian.mat');
%vpa(x,289);
%vpa(weight,289);

% fid = fopen('F:\Git\Cetim\ep_a_04\Acqui_CV_1.txt');
fid = fopen('/home/ma/MATLAB/Cetim/ep_a_04/Acqui_CV_1.txt');
area=2.27*9.95*1e-6; %meter square
NF=2500000;
run('Referece_scalar_iteration');
xlwrite('a_fitting.xls',nF,1,'D04');
xlwrite('a_fitting.xls',(NF-nF)/NF,1,'E04');
