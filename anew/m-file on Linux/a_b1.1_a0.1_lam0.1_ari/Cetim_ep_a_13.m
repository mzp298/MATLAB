clear;clc;
close all
format long e
run('Parameters_ep_a_b')
load('gaussian.mat');
%vpa(x,289);
%vpa(weight,289);

% fid = fopen('F:\Git\Cetim\ep_a_13\Acqui_CV.txt');
fid = fopen('/home/ma/MATLAB/Cetim/ep_a_13/Acqui_CV.txt');
area=2.28*9.99*1e-6; %meter square
NF=10763158;
run('Referece_scalar_iteration');
xlwrite('a_fitting.xls',nF,1,'D13');
xlwrite('a_fitting.xls',(NF-nF)/NF,1,'E13');
