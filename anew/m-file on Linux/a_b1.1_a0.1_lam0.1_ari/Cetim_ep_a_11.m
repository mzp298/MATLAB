clear;clc;
close all
format long e
run('Parameters_ep_a_b')
load('gaussian.mat');
%vpa(x,289);
%vpa(weight,289);

% fid = fopen('F:\Git\Cetim\ep_a_11\Acqui_CV.txt');
fid = fopen('/home/ma/MATLAB/Cetim/ep_a_11/Acqui_CV.txt');
area=2.29*10.01*1e-6; %meter square
NF=10131579;
run('Referece_scalar_iteration');
xlwrite('a_fitting.xls',nF,1,'D11');
xlwrite('a_fitting.xls',(NF-nF)/NF,1,'E11');
