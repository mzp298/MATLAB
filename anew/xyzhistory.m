clear;clc;close all;

figure(2)
load('FX_RAVG.mat')
x=plot(signal.data,'b');
hold on
load('FY_RAVG.mat')
y=plot(signal.data,'r')
load('FZ_RAVG.mat')
z=plot(signal.data,'m')

grid on;
grid minor;
 hLegend=legend({'rear to front(X)','left to right(Y)','top to bottom(Z)'},'FontSize',25,'FontWeight','bold');
 set([hLegend, gca], 'FontSize', 20)
set(hLegend,'Box','on');
set(hLegend,'EdgeColor',[1 1 1]); %set the edge colour of the legend to white 
% hTitle =title('Force history on a car suspension arm in 3 directions in cartesian coordinate');
% hTitle =title('Stress history of top to botom direction(Z direction)','Fontsize' ,25);
hXLabel =xlabel('number of recorded points with frenquency f=256(s^{-1})' ,'Fontsize' ,25,'FontWeight' , 'bold');
 hYLabel =ylabel('Force(N)', 'Fontsize' ,25,'FontWeight' , 'bold');
 
 set(gca, 'FontName', 'Helvetica')
set([hXLabel, hYLabel], 'FontName', 'AvantGarde')
set([hLegend, gca], 'FontSize', 20)
set([hXLabel, hYLabel], 'FontSize', 20)


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
set(gcf, 'PaperPosition', [0 0 1080 600]); %set(gcf,'PaperPosition',[left,bottom,width,height])
 saveas(gcf,'F:\Git\Anew\figures\xyz_suspension.png');
  saveas(gcf,'F:\Git\PhDreport\BEAMER\UTC2017\figures\xyz_suspension.png');