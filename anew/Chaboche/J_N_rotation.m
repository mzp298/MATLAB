clear;clc;close all;
data=xlsread('f:\stress values rotating cantilever');
A1=data(5,4);
B1=data(5,5);
lam=data(3,1);
miu=data(3,2);
tt=data(7,1);
ff=data(7,2);
ac=(tt-ff/sqrt(3))/(ff/3);
bc=tt;
ad=(3*tt)/(ff)-(1.5);
bd=tt;
sigu=1.67e6;
r=0.1;

sig1=[(lam*(4*A1*(r^2)+2*B1)+2*miu*(3*A1*(r^2)+B1)) 0 0;0 (lam*(4*A1*...
(r^2)+2*B1)+2*miu*(A1*(r^2)+B1)) 0;0 0 (lam*(4*A1*(r^2)+2*B1));];
p1=1/3*sum(diag(sig1));
%p is independent of t, so pmax=p
S1=2*sig1-(1/3*sum(diag(sig1)))*diag([1,1,1]);
sqrj1=1/2*sqrt(1/2*(S1(1,1)^2+S1(2,2)^2+S1(3,3)^2+2*(S1(1,2)^2)+...
2*(S1(1,3)^2)+2*(S1(2,3)^2)));
pm1=p1;
%---------------Modification with Crossland and DangVan-------------
cross1=sqrj1+ac*pm1-bc;
tau1=1/2*(sig1(2,2)-sig1(3,3));
dangvan1=tau1+ad*pm1-bd;
beta=6;
M=ff*(1-3*pm1/sigu);
a=0.1;
b=5*ac/(3*ff);%we set parameter '3*b*ff' 5 times as 'ac' in Crossland

alphas=1-a*(sqrj1-ff*(1-3*b*pm1))/(sigu-2*sqrj1);
alphac=1-a*(cross1)/(sigu-2*sqrj1);
alphad=1-a*(dangvan1)/(sigu-2*sqrj1);

NF=0:1e5:1e9;
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
hTitle = title('A_{||}-N curve in ROTATION using Chaboche model with different criteria' ,'Fontsize' ,30);
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


% saveas(gcf,'F:\Git\PhDreport\5thesis\figures\JNrotation.png');

sp=actxserver('SAPI.SpVoice');
sp.Speak('I have finished the job you gave me and thank you for that');
