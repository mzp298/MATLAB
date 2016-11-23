clear;
clc

s=1:0.001:30;
beta=1.5;
c=beta-1;
p=c*s.^-beta;
cum=1-s.^(1-beta);
hold on
%loglog(s,p,'r')
plot(s,p,'LineStyle', 'none','LineWidth', 1,'Marker', 'o', 'MarkerSize', 1, ...
      'MarkerEdgeColor','none', 'MarkerFaceColor',[238 180 34]/255);
h=area(s,p,'LineStyle', 'none','LineWidth', 1);
set(h,'FaceColor',[238 180 34]/255);
set(gca,'FontSize',25);
set(gca,'xtick',1:5:30,'xticklabel',1:5:30) 
%xtick shows the ticks, xticklabel gives the tick names.
set(gca,'xlim',[0 30]);
set(gca,'ylim',[0 0.6]);
grid on;
%grid minor;
hTitle = title('Weakening scales distribution curve','Fontsize',35);
hXLabel = xlabel('Scale s','Fontsize',30);
hYLabel = ylabel('Probability distribution function P(s)','Fontsize',30);
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
set(gcf, 'PaperPosition', [0 0 1920 1080]); 

% saveas(gcf,'F:\Git\PhDreport\4Anew\figures\ps.png');
% saveas(gcf,'F:\Git\Anew\figures\ps.png');