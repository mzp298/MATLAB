clc;
close all;

s=1:0.001:20;
beta=5;
c=beta-1;
p=c*s.^-beta;
cum=1-s.^(1-beta);
hold on
%loglog(s,p,'r')
plot(s,p,'LineStyle', 'none','LineWidth', 1,'Marker', 'o', 'MarkerSize', 1, ...
      'MarkerEdgeColor','none', 'MarkerFaceColor',[138 43 226]/255);
h=area(s,p,'LineStyle', 'none','LineWidth', 1);
set(h,'FaceColor',[138 43 226]/255);
set(gca,'FontSize',25);
set(gca,'xtick',1:5:30,'xticklabel',1:5:30) 
%xtick shows the ticks, xticklabel gives the tick names.
set(gca,'xlim',[0 max(s)]);
set(gca,'ylim',[0 max(p)]);
grid on;
%grid minor;
hXLabel = xlabel('Scale s','Fontsize',35);
hYLabel = ylabel('Probability distribution function P(s)','Fontsize',35);
set(gca, 'FontName', 'Helvetica')
set([hXLabel, hYLabel], 'FontName', 'AvantGarde')
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
set(gcf, 'PaperPosition', [0 0 1920 1080]); 

% saveas(gcf,'F:\Git\PhDreport\4Anew\figures\ps2.png');
saveas(gcf,'F:\Git\Anew\figures\ps2.png');