clear;clc;close all;
D=0:0.01:1;
gam=2;
Dhat=1-(1-D).^(gam+1);
figure(1);
Dhat=plot(D,Dhat,'-k','LineWidth',4);
axis equal;
axis([0 1 0 1]);
grid on;

hXLabel = xlabel('D');
hYLabel =ylabel({'$\hat{D}$'},'Interpreter','latex');
set(gca,'XTick',0:0.1:1);
set(gca,'YTick',0:0.1:1);

% Adjust font
set(gca, 'FontName', 'Helvetica')
set([ hXLabel, hYLabel], 'FontName', 'AvantGarde')
set([hXLabel, hYLabel, gca], 'FontSize',30)

% Adjust axes properties
set(gcf,'color','w'); %set figure background transparent
set(gca,'color','w'); %set axis transparent
% Maximize print figure
set(gcf,'outerposition',get(0,'screensize'));
set(gcf, 'PaperPositionMode', 'manual');
set(gcf, 'PaperUnits', 'points'); %[ {inches} | centimeters | normalized | points ]
set(gcf, 'PaperPosition', [0 0 1000 1000]); %set(gcf,'PaperPosition',[left,bottom,width,height])
saveas(gcf,'F:\Git\Anew\figures\Dhat.png');

