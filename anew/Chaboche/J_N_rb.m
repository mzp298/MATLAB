clear;clc;close all;
data=xlsread('f:\stress values rotating cantilever');
A=-1.19e-7;
B=3.34e-7;
lam=data(3,1);
miu=data(3,2);
tt=300e6;
ff=400e6;
ac=(tt-ff/sqrt(3))/(ff/3);
ad=(3*tt)/(ff)-(1.5);
bc=tt;
bd=tt;
w=5; %---------------to be adjusted-----------------------
F1=6e8; %---------------to be adjusted-----------------------

L=50;
I=127.2345025;
sigu=1000e6;
%r=0.5;
z=0 ; 
r=3 ;

syms t ;
sig1=[(lam*(4*A*(r^2)+2*B)+2*miu*(3*A*(r^2)+B)) 0 0;0 (lam*(4*A*...
(r^2)+2*B)+2*miu*(A*(r^2)+B)) 0;0 0 (lam*(4*A*(r^2)+2*B))-F1*(L-...
z)/I*r*sin(w*t+pi/2);];

p1=1/3*sum(diag(sig1));

S1=sig1-(1/3*sum(diag(sig1)))*diag([1,1,1]);
sqrj1=1/2*sqrt(1/2*(S1(1,1)^2+S1(2,2)^2+S1(3,3)^2+2*(S1(1,2)^2)+...
2*(S1(1,3)^2)+2*(S1(2,3)^2)));

b=5*ac/(3*ff);%we set parameter '3*b*ff' 5 times as 'ac' in Crossland

%---------------Modification with Crossland and DangVan-------------
t=0:1e-2:3;

p1=subs(p1);
pm=vpa((max(p1)+min(p1))/2,5);
sqrj1=subs(sqrj1);
sqrj1=vpa(max(sqrj1),5);

pM=vpa(max(p1),5);
cross1=sqrj1+ac*pM-bc;

tau1 = 1/2*(sig1(3,3)-sig1(1,1));
dangvan1=tau1+ad*p1-bd;
dangvan1=subs(dangvan1);
dangvan1=vpa(max(dangvan1),5);

beta=6;
M=ff*(1-3*pm/sigu);
a=0.1;
b=1e-8*ac/(3*ff);%we set parameter '3*b*ff' 5 times as 'ac' in Crossland

alphas=1-a*(sqrj1-ff*(1-3*b*pm))/(sigu-2*sqrj1);
alphac=1-a*(cross1)/(sigu-2*sqrj1);
alphad=1-a*(dangvan1)/(sigu-2*sqrj1);

NF=1:1e3:1e6;
As=(NF*(beta+1)*(1-alphas)).^(-1/beta)*M;
Ac=(NF*(beta+1)*(1-alphac)).^(-1/beta)*M;
Ad=(NF*(beta+1)*(1-alphad)).^(-1/beta)*M;


PAs=loglog(NF,As,'g','LineStyle', ':','LineWidth', 7);
hold on
PAc=loglog(NF,Ac,'r','LineStyle', '-','LineWidth', 7);
PAd=loglog(NF,Ad,'b','LineStyle', '--','LineWidth', 7);

grid on;
grid minor;
set(gca ,'FontSize',30);
hXLabel = xlabel('N_F' ,'Fontsize' ,30);
hYLabel = ylabel('A_{||}', 'Fontsize' ,30);
hTitle = title('A_{||}-N curve in ROTATIVE BENDING using Chaboche model with different criteria' ,'Fontsize' ,30);
set(hTitle, 'FontSize', 30, 'FontWeight' , 'bold')
hLegend=legend('Sines criteria','Crossland criteria','Dangvan criteria',...
    'Location','best');
set(hLegend, 'FontSize', 20)
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
set(gcf, 'PaperPosition', [0 0 1280 800]); %set(gcf,'PaperPosition',[left,bottom,width,height])


saveas(gcf,'F:\Git\PhDreport\5thesis\figures\JNRB.png');

sp=actxserver('SAPI.SpVoice');
sp.Speak('Dont forget it');

