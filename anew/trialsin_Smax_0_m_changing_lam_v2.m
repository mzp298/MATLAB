% Program to get the Gauss-Legendre Quadrature results (Vectorized)
clear;clc;close all;
format long e
load('gaussian.mat');

E=191e9;               %Young's modulus
nu=0.38;                 %poisson's ratio
k=1e9;                  %hardening parameter
b=1.5;                      %weakening scales distribution exponent (between 1 and 2)
fb=1.1;
a=0.1;
delta_alp=1e-3;
y=1080e6;            %macroscopic yield stress
load=7e8;            %cyclic load
loadtensor= [load 0 0;0 0 0;0 0 0];
stepnumber=100;        %devide one cycle in 200 parts
cycles=2;
delta_alp=1e-4;
f=50;                            %frequency of load
steptime=1/f/stepnumber;
W0=5e8;             %dissipated energy to failure per unit volume
lamratio=1;
lamplus=0.6;
lamminus=lamratio.*lamplus;
m=3e8;                   % mean stress
% m=0;                   % mean stress
sigm=m;
scentre=[2*sigm/3            0                0 ;...
    0             -sigm/3                0 ;...
    0                       0      -sigm/3];
hydrofix=1/3*(sum(diag(loadtensor)));
sig=loadtensor-hydrofix*eye(3); %mean stress does not change deviatoric stress!!!!!!!!

%---------------------3 Numerical method-----------------------------
D=1e-16;
n=1;
tensor = [m+load*sind(n*360/stepnumber) 0 0 ;...
    0 0 0 ;...
    0 0 0 ];
run('Damiter1.m')
D=D+D^alp*W/W0;
tic;
while n<2*stepnumber
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
    D=D+D^alp(n+1)*W/W0;
    figure(1);
    hold on;
    yield1=plot (n+1,yield(n+1)*s(3).^-1, 'LineStyle', 'none','LineWidth', 1, 'Marker', 'o', 'MarkerSize', 8, ...
        'MarkerEdgeColor',  'none', 'MarkerFaceColor' , 'c');
    yield1n= plot (n+1,-yield(n+1)*s(3).^-1, 'LineStyle', 'none','LineWidth', 1, 'Marker', 'o', 'MarkerSize', 8, ...
        'MarkerEdgeColor',  'none', 'MarkerFaceColor' , 'c');
    Trial1=plot (n+1,sign(trial11(3))*Smaxtrial(3),'LineStyle', 'none','LineWidth', 1,'Marker', '^', 'MarkerSize', 10, ...
        'MarkerEdgeColor','r', 'MarkerFaceColor','r');
    Sb1=plot (n+1,sign(Sb11(3))*normSb(3),'LineStyle', 'none','LineWidth', 1,'Marker', 'v', 'MarkerSize', 10, ...
        'MarkerEdgeColor','g', 'MarkerFaceColor','g');
    sig=plot (n+1,tensor(1,1),'LineStyle', 'none','LineWidth', 1,'Marker', 's', 'MarkerSize', 11, ...
        'MarkerEdgeColor','none', 'MarkerFaceColor',[238 18 137]/255);
    yield8=plot (n+1,yield(n+1)*s(10).^-1,'LineStyle', 'none','LineWidth', 1,'Marker', 'o', 'MarkerSize', 8, ...
        'MarkerEdgeColor', 'none', 'MarkerFaceColor', 'b');
    yield8n=plot (n+1,-yield(n+1)*s(10).^-1,'LineStyle', 'none','LineWidth', 1,'Marker', 'o', 'MarkerSize', 8, ...
        'MarkerEdgeColor', 'none', 'MarkerFaceColor', 'b');
    Trial8=plot (n+1,sign(trial11(10))*Smaxtrial(10),'LineStyle', 'none','LineWidth', 1,'Marker', '^', 'MarkerSize', 10, ...
        'MarkerEdgeColor', [1 0.5 0], 'MarkerFaceColor',[1 0.5 0]);
    Sb8=plot (n+1,sign(Sb11(10))*normSb(10),'LineStyle', 'none','LineWidth', 1,'Marker', 'v', 'MarkerSize', 10, ...
        'MarkerEdgeColor','k', 'MarkerFaceColor','k');
    
    n=n+1;
end
toc;
n

Nf=n*stepnumber^-1;
disp(['Cycles to failure is ' num2str(Nf) ' cycles.']);

% %---------------------in loop 2 scales plot settings-----------------------------
figure(1);
grid on;
grid minor;
axis([0 n -10e8 10e8]);
set(gca ,'FontSize',30);
hXLabel = xlabel('Time step(100 in unit cycle)' ,'Fontsize' ,30);
hYLabel = ylabel('Stress(Pa)', 'Fontsize' ,30);
hTitle = title('Microscopic stress evolution at 2 scales' ,'Fontsize' ,30);
set(hTitle, 'FontSize', 30, 'FontWeight' , 'bold')
hLegend=legend([sig,yield1,Sb1,Trial1,yield8,Sb8,Trial8],'{$\underline{\underline{\Sigma}}_{bending}$}',...
    '{$(\sigma_y-\lambda\Sigma_H)/s_{3} \quad\;\,\,  at \; scale\; s_{3}$}',...
    '{$||S-b||   \quad\quad\qquad\;\;\,  at \; scale\; s_{3}$}',...
    '{$||S-b||_{trial}   \quad\qquad\,   at \; scale\; s_{3}$}',...
    '{$(\sigma_y-\lambda\Sigma_H)/s_{10} \quad\;\,   at \; scale\; s_{10}$}',...
    '{$||S-b||   \quad\quad\qquad\;\;\,   at \; scale\; s_{10}$}',...
    '{$||S-b||_{trial}   \quad\qquad\,   at \; scale\; s_{10}$}',...
    'Location','northeast');
set(hLegend, 'Interpreter', 'latex')
set(hLegend, 'FontSize',20)
set(hLegend,'Box','on');
set(hLegend,'EdgeColor',[1 1 1]); %set the edge colour of the legend to white
% Adjust font
set(gca, 'FontName', 'Helvetica')
set([hTitle, hXLabel, hYLabel], 'FontName', 'AvantGarde')
% Adjust axes properties
set(gca, 'Box', 'off', 'TickDir', 'out', 'TickLength', [.02 .02], ...
    'XMinorTick', 'on', 'YMinorTick', 'on', 'YGrid', 'on', ...
    'XColor', [.3 .3 .3], 'YColor', [.3 .3 .3], ...
    'LineWidth', 1)
set(gcf,'color','w'); %set figure background transparent
set(gca,'color','w'); %set axis transparent
% Maximize print figure
set(gcf,'outerposition',get(0,'screensize'));
set(gcf, 'PaperPositionMode', 'manual');
set(gcf, 'PaperUnits', 'points'); %[ {inches} | centimeters | normalized | points ]
set(gcf, 'PaperPosition', [0 0 1280 720]); %set(gcf,'PaperPosition',[left,bottom,width,height])
% saveas(gcf,'F:\Git\Anew\figures\trialsin_m.png');
%  saveas(gcf,'F:\Git\Anew\figures\trialsin_0.png');
%
sp=actxserver('SAPI.SpVoice');
sp.Speak('Patrick Le tallec');


