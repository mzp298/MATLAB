% integration of  N over D
clc;
clear;
a=0.5;
b=3.5;
A= 7.7606E5;
M=1e8;
D=0:0.01:1;

N=(1 - (1 - D).^(b + 1)).^(1-a)*(A/M).^-b/(1+b)/(1-a);
dn=plot(N,D,'LineWidth', 10);
hXLabel =xlabel('Number of cycles N');
hYLabel=ylabel('Damage D');


grid on;
grid minor;
% hTitle =title('S-N curve at constant amplitude stress with R=-1 in our model with different $W_F$');
hXLabel = xlabel('N_F');
hYLabel =ylabel('Damage')

% Adjust font
set(gca, 'FontName', 'Helvetica','FontSize',40)
set([hXLabel, hYLabel], 'FontName', 'AvantGarde')
set([hXLabel, hYLabel], 'FontSize',40)

% Adjust axes properties
set(gcf,'color','w'); %set figure background transparent
set(gca,'color','w'); %set axis transparent
% Maximize print figure
set(gcf,'outerposition',get(0,'screensize'));
set(gcf, 'PaperPositionMode', 'manual');
set(gcf, 'PaperUnits', 'points'); %[ {inches} | centimeters | normalized | points ]
set(gcf, 'PaperPosition', [0 0 1200 1000]); %set(gcf,'PaperPosition',[left,bottom,width,height])
%  saveas(gcf,'F:\Git\PhDreport\5thesis\figures\D-N.png');