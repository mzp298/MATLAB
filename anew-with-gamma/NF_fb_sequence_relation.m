% Program to get the Gauss-Legendre Quadrature results (Vectorized)
clear;clc;

dbstop if error
format long e

[x]= [0.999305042	0.996340117	0.991013371	0.983336254	0.973326828	0.9610088	0.946411375	0.929569172	0.910522137...
    0.889315446	0.865999398	0.840629296	0.813265315	0.783972359	0.752819907	0.71988185	0.685236313	0.648965471...
    0.611155355	0.571895646	0.531279464	0.489403146	0.446366017	0.402270158	0.357220158	0.311322872	0.264687162...
    0.217423644	0.16964442	0.121462819	0.072993122	0.024350293	-0.024350293	-0.072993122	-0.121462819	-0.16964442...
    -0.217423644	-0.264687162	-0.311322872	-0.357220158	-0.402270158	-0.446366017	-0.489403146	-0.531279464...
    -0.571895646	-0.611155355	-0.648965471	-0.685236313	-0.71988185	-0.752819907	-0.783972359	-0.813265315...
    -0.840629296	-0.865999398	-0.889315446	-0.910522137	-0.929569172	-0.946411375	-0.9610088	-0.973326828...
    -0.983336254	-0.991013371	-0.996340117	-0.999305042];
[weight]=[0.001783281	0.004147033	0.006504458	0.00884676	0.011168139	0.013463048	0.01572603	0.017951716	0.020134823...
    0.022270174	0.024352703	0.02637747	0.028339673	0.030234657	0.032057928	0.033805162	0.035472213	0.037055129	0.038550153...
    0.039953741	0.041262563	0.042473515	0.043583725	0.044590558	0.045491628	0.046284797	0.046968183	0.047540166	0.047999389...
    0.048344762	0.048575467	0.048690957	0.048690957	0.048575467	0.048344762	0.047999389	0.047540166	0.046968183	0.046284797...
    0.045491628	0.044590558	0.043583725	0.042473515	0.041262563	0.039953741	0.038550153	0.037055129	0.035472213	0.033805162...
    0.032057928	0.030234657	0.028339673	0.02637747	0.024352703	0.022270174	0.020134823	0.017951716	0.01572603	0.013463048...
    0.011168139	0.00884676	0.006504458	0.004147033	0.001783281];
% [x]=xlsread('Gauss-Legendre Quadrature','Sheet1','b1:z1');
% [weight]=xlsread('Gauss-Legendre Quadrature','Sheet1','b2:z2');

E=191e9;               %Young's modulus
nu=0.38;                 %poisson's ratio
k=1e9;                  %hardening parameter
b=3;                      %weakening scales distribution exponent (between 1 and 2)
y=1080e6;            %macroscopic yield stress
sigu=1200e6;             %ultimite stress
ff=690e6;              %bending fatigue limit
tt=428e6;                  %torsion fatigue limit

gam=0.1;              %material parameter from Chaboche law(Wohler curve exponent)

load=7e8;            %cyclic load
loadtensor= [load 0 0;0 0 0;0 0 0];
stepnumber=32;        %devide one cycle in 200 parts
f=50;                            %frequency of load
steptime=1/f/stepnumber;
delta=(b+1)/(b-1);
alp=0.5;
WF=5e8;             %dissipated energy to failure per unit volume
n0=1;                   %number of initial local defects
lam=0.3;               %hydrostatic pressure sensitivity
m=3e8;                   % mean stress
hydrofix=1/3*(sum(diag(loadtensor))); 
dev=loadtensor-hydrofix*eye(3); %mean stress does not change deviatoric stress!!!!!!!!
a=0.1;                %sensitivity of sequence effect(control alp>0)
X=1.065;                 %Major damage threshold (meaning S_{max} must be greater than X/2*yield to activate magnification)
pb=b;                %Major damage power

%---------------------2 Cyclic load calculation-----------------------------
Dcyc(1)=0;
n=1;
Gcyc = (1 - (1 - Dcyc(1)).^(gam + 1)).^(1-alp);
Smax=norm(dev,'fro'); 
yield=y-lam*1/3*m; %mean value of hydro is 1/3*m
Wcyc=4*(E-k)*(1+nu)*(b-1)/(E*(E+k*nu)*b*(b+1))*Smax.^(b+1)*yield.^(1-b);
sequence=((Smax*yield^-1)*(X-Smax*yield^-1)^-1)^pb;
NF1=1/((1+gam)*a*sequence)
NF2=WF/n0/Wcyc
NF=NF1*NF2
%  sequence<1 b=1.5 NF1 =  2.0
%                              NF2 = 1.2
%                              NF =  26
%                      b=3   NF1 = 4.7
%                               NF2 =4.1
%                               NF =  190
% sequence>1 b=1.5  NF1 =8.5
%                              NF2 =5.4
%                              NF =46
%                      b=3  NF1 =7.9 %decrease!!!!!!!!!!!!
%                              NF2 =10
%                              NF =85