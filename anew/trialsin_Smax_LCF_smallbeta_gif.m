% Program to get the Gauss-Legendre Quadrature results (Vectorized)
clear;clc;close all;
format long e
load('gaussian.mat');

E=191e9;               %Young's modulus
nu=0.38;                 %poisson's ratio
k=1e9;                  %hardening parameter
b=2;                      %weakening scales distribution exponent (between 1 and 2)
fb=1.1;
a=0.1;
delta_alp=1e-3;
y=1080e6;            %macroscopic yield stress
load=8.5e8;            %cyclic load
loadtensor= [load 0 0;0 0 0;0 0 0];
stepnumber=100;        %devide one cycle in 200 parts
cycles=2;
delta_alp=1e-4;
f=50;                            %frequency of load
steptime=1/f/stepnumber;
W0=5e8;             %dissipated energy to failure per unit volume
lamplus=0.7;
lamminus=0.2;
% m=3e8;                   % mean stress
m=0;                   % mean stress
sigm=m;
scentre=[2*sigm/3            0                0 ;...
                0              -sigm/3                0 ;...
                0                         0      -sigm/3];
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
figure('name','LCF');
filename = 'trialsin_Smax_LCF.gif';
grid on;
grid minor;
axis([0 2*stepnumber -1.1*y  1.1*y]);
set(gca ,'FontSize',30);

% Adjust axes properties
set(gca, 'Box', 'off', 'TickDir', 'out', 'TickLength', [.02 .02], ...
    'XMinorTick', 'on', 'YMinorTick', 'on', 'YGrid', 'on', ...
    'XColor', [.3 .3 .3], 'YColor', [.3 .3 .3], ...
    'LineWidth', 1);
set(gcf,'color','w'); %set figure background transparent
set(gca,'color','w'); %set axis transparent
% Maximize print figure
set(gcf,'outerposition',get(0,'screensize'));
set(gcf, 'PaperPositionMode', 'manual');
set(gcf, 'PaperUnits', 'points'); %[ {inches} | centimeters | normalized | points ]
set(gcf, 'PaperPosition', [0 0 1280 720]); %set(gcf,'PaperPosition',[left,bottom,width,height])
hXLabel = xlabel('Time step(100 in unit cycle)' ,'Fontsize' ,30);
hYLabel = ylabel('Stress(Pa)', 'Fontsize' ,30);
hTitle = title('Stress intensity evolution at different scales for LCF' ,'Fontsize' ,30);
set(hTitle, 'FontSize', 30, 'FontWeight' , 'bold');

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
    
    hold on;
    yield0=plot (n+1,yield(n+1), 'LineStyle', 'none','LineWidth', 1, 'Marker', 'o', 'MarkerSize', 8, ...
        'MarkerEdgeColor',  'none', 'MarkerFaceColor' , [0 0 139]/255);
    yield0n= plot (n+1,-yield(n+1), 'LineStyle', 'none','LineWidth', 1, 'Marker', 'o', 'MarkerSize', 8, ...
        'MarkerEdgeColor',  'none', 'MarkerFaceColor' , [0 0 139]/255);
    sig=plot (n+1,tensor(1,1),'LineStyle', 'none','LineWidth', 1,'Marker', 's', 'MarkerSize', 11, ...
        'MarkerEdgeColor','none', 'MarkerFaceColor',[238 18 137]/255);
    
    yield1=plot (n+1,yield(n+1)*s(900).^-1, 'LineStyle', 'none','LineWidth', 1, 'Marker', 'o', 'MarkerSize', 8, ...
        'MarkerEdgeColor',  'none', 'MarkerFaceColor' , 'c');
    yield1n= plot (n+1,-yield(n+1)*s(900).^-1, 'LineStyle', 'none','LineWidth', 1, 'Marker', 'o', 'MarkerSize', 8, ...
        'MarkerEdgeColor',  'none', 'MarkerFaceColor' , 'c');
    Trial1=plot (n+1,sign(trial11(900))*Smaxtrial(900),'LineStyle', 'none','LineWidth', 1,'Marker', '^', 'MarkerSize', 10, ...
        'MarkerEdgeColor','r', 'MarkerFaceColor','r');
    Sb1=plot (n+1,sign(Sb11(900))*normSb(900),'LineStyle', 'none','LineWidth', 1,'Marker', 'v', 'MarkerSize', 10, ...
        'MarkerEdgeColor','g', 'MarkerFaceColor','g');
    
    yield8=plot (n+1,yield(n+1)*s(700).^-1,'LineStyle', 'none','LineWidth', 1,'Marker', 'o', 'MarkerSize', 8, ...
        'MarkerEdgeColor', 'none', 'MarkerFaceColor', [28 134 238]/255);
    yield8n=plot (n+1,-yield(n+1)*s(700).^-1,'LineStyle', 'none','LineWidth', 1,'Marker', 'o', 'MarkerSize', 8, ...
        'MarkerEdgeColor', 'none', 'MarkerFaceColor', [28 134 238]/255);
    Trial8=plot (n+1,sign(trial11(700))*Smaxtrial(700),'LineStyle', 'none','LineWidth', 1,'Marker', '^', 'MarkerSize', 10, ...
        'MarkerEdgeColor', [1 0.5 0], 'MarkerFaceColor',[1 0.5 0]);
    Sb8=plot (n+1,sign(Sb11(700))*normSb(700),'LineStyle', 'none','LineWidth', 1,'Marker', 'v', 'MarkerSize', 10, ...
        'MarkerEdgeColor','k', 'MarkerFaceColor','k');
    
    yield9=plot (n+1,yield(n+1)*s(550).^-1,'LineStyle', 'none','LineWidth', 1,'Marker', 'o', 'MarkerSize', 8, ...
        'MarkerEdgeColor', 'none', 'MarkerFaceColor', [142 229 238]/255);
    yield9n=plot (n+1,-yield(n+1)*s(550).^-1,'LineStyle', 'none','LineWidth', 1,'Marker', 'o', 'MarkerSize', 8, ...
        'MarkerEdgeColor', 'none', 'MarkerFaceColor', [142 229 238]/255);
    Trial9=plot (n+1,sign(trial11(550))*Smaxtrial(550),'LineStyle', 'none','LineWidth', 1,'Marker', '^', 'MarkerSize', 10, ...
        'MarkerEdgeColor', 'none', 'MarkerFaceColor',[238 99 99]/255);	
    Sb9=plot (n+1,sign(Sb11(550))*normSb(550),'LineStyle', 'none','LineWidth', 1,'Marker', 'v', 'MarkerSize', 10, ...
        'MarkerEdgeColor','none', 'MarkerFaceColor',[154 205 50]/255);
    
    % %---------------------in loop 2 scales plot settings-----------------------------
    hLegend=legend([yield0,sig,yield1,Sb1,Trial1,yield8,Sb8,Trial8,yield9,Sb9,Trial9],'{$\sigma_y$}',...
        '{$\underline{\underline{\Sigma}}_{bending}$}',...
        '{$(\sigma_y-\lambda\Sigma_H)/s_{3} \quad\;\,\,  at \; scale\; s_{3}$}',...
        '{$||S-b||   \quad\quad\qquad\;\;\,  at \; scale\; s_{3}$}',...
        '{$||S-b||_{trial}   \quad\qquad\,   at \; scale\; s_{3}$}',...
        '{$(\sigma_y-\lambda\Sigma_H)/s_{10} \quad\;\,   at \; scale\; s_{10}$}',...
        '{$||S-b||   \quad\quad\qquad\;\;\,   at \; scale\; s_{10}$}',...
        '{$||S-b||_{trial}   \quad\qquad\,   at \; scale\; s_{10}$}',...
        '{$(\sigma_y-\lambda\Sigma_H)/s_{20} \quad\;\,   at \; scale\; s_{20}$}',...
        '{$||S-b||   \quad\quad\qquad\;\;\,   at \; scale\; s_{20}$}',...
        '{$||S-b||_{trial}   \quad\qquad\,   at \; scale\; s_{20}$}',...
        'Location','northeast');
    set(hLegend, 'Interpreter', 'latex');
    set(hLegend, 'FontSize',20);
    set(hLegend,'Box','on');
    set(hLegend,'EdgeColor',[1 1 1]); %set the edge colour of the legend to white
    % Adjust font
    set(gca, 'FontName', 'Helvetica');
    set([hTitle, hXLabel, hYLabel], 'FontName', 'AvantGarde');
    
    %%---------------to generate gif slide by slide---------------------
    frame = getframe(1);
    im = frame2im(frame);
    [A,map] = rgb2ind(im,256);
    if n == 1;
        imwrite(A,map,filename,'gif','LoopCount',Inf,'DelayTime',0.001);
    else
        imwrite(A,map,filename,'gif','WriteMode','append','DelayTime',0.001);
    end
    
    n=n+1;
end
toc;

saveas(gcf,'trialsin_Smax_LCF.png');
sp=actxserver('SAPI.SpVoice');
sp.Speak('Patrick Le tallec');


