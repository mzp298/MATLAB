%----------change parameters to colaborate with CETIM tests-------------
clear;clc;close all;


b=1.07;                    %weakening scales distribution exponent
y=230e6;           %macroscopic yield stress
lam=0.1;               %hydrostatic pressure sensitivity
a=0.4;
W0=5.4e8;            %dissipated energy to failure per unit volume
E=72e9;              %Young's modulus
fb=b;                %Major damage power
k=6e8;                 %hardening parameter
nu=0.3;                     %poisson's ratio
ari=100;  %divide a reverse into ''ari'' parts
delta_alp=0.01;
repetition=26316;
%save('CETIM.mat','a','b','lam','W0','E','fb','k','nu','y','ari','delta_alp','repetition') 
%% Initialisation of POI Libs
% % Add Java POI Libs to matlab javapath
javaaddpath('poi_library/poi-3.8-20120326.jar');
javaaddpath('poi_library/poi-ooxml-3.8-20120326.jar');
javaaddpath('poi_library/poi-ooxml-schemas-3.8-20120326.jar');
javaaddpath('poi_library/xmlbeans-2.3.0.jar');
javaaddpath('poi_library/dom4j-1.6.1.jar');
javaaddpath('poi_library/stax-api-1.0.1.jar');


