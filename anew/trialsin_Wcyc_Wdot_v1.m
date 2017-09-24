% Program to get the Gauss-Legendre Quadrature results (Vectorized)
clear;clc;close all;
format long e

load('gaussian.mat');

E=191e9;               %Young's modulus
nu=0.38;                 %poisson's ratio
k=1e9;                  %hardening parameter
b=1.5;                      %weakening scales distribution exponent (between 1 and 2)
y=1080e6;            %macroscopic yield stress
fb=1.1;
load=6e8;            %cyclic load
stepnumber=100;        %devide one cycle in 200 parts
a=0.5;
delta_alp=1e-3;
W0=1.5e6;             %dissipated energy to failure per unit volume
lamratio=1;
lamplus=0.3;
lamminus=lamratio.*lamplus;
m=0;                   % mean stress
% m=3e8;                   % mean stress
sigm=m;
scentre=[2*sigm/3            0                0 ;...
                 0             -sigm/3                0 ;...
                 0                       0      -sigm/3];
%---------------------3 Numerical method with alp(Smin)-----------------------------
D2= ones(1,1e6)*1e-16; %Pre-allocate memory for vectors
n=1;       %initial recording point
tensor = [m+load*sind(n*360/stepnumber) 0 0 ;...
    0 0 0 ;...
    0 0 0 ];
run('Damiter1.m')
D2(1)=D2(1)+D2(1)^alp(1)*W/W0;
tic;
% while n<2*stepnumber
while D2(n)<1
    tensor = [m+load*sind(n*360/stepnumber) 0 0 ;...
        0 0 0 ;...
        0 0 0 ];
    hydro(n)=1/3*trace(tensor);
    dev1=tensor-hydro(n)*eye(3)-scentre;
    dev11=dev1(1,1); dev12=dev1(1,2); dev13=dev1(1,3);
    dev21=dev1(2,1); dev22=dev1(2,2); dev23=dev1(2,3);
    dev31=dev1(3,1); dev32=dev1(3,2); dev33=dev1(3,3);
    tensor = [m+load*sind((n+1)*360/stepnumber) 0 0 ;...
        0 0 0 ;...
        0 0 0 ];
    run('Damiter2.m')
%     figure(1);
%     hold on;
%     W_change_alp=plot ((n+1),W,'LineStyle', 'none','LineWidth', 1, 'Marker', 's', 'MarkerSize',8, ...
%         'MarkerEdgeColor',  [255 193 37]/255, 'MarkerFaceColor' , [255 193 37]/255);
%     sig_scale=plot (n+1,tensor(1,1)/4.2e5,'LineStyle', 'none','LineWidth', 1,'Marker', 's', 'MarkerSize', 11, ...
%         'MarkerEdgeColor','none', 'MarkerFaceColor',[238 18 137]/255);
    
    D2(n+1)=D2(n)+D2(n)^alp(n+1)*W/W0;
    n=n+1;
end

toc;
n

disp(['Mean alpha= ' num2str(mean(alp)) ' .']);


sp=actxserver('SAPI.SpVoice');
sp.Speak('finished');
