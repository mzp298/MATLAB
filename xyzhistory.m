load('FX_RAVG.mat')
x=plot(signal.data,'b')
hold on
load('FY_RAVG.mat')
y=plot(signal.data,'r')
load('FZ_RAVG.mat')
z=plot(signal.data,'m')

grid on;
grid minor;
 hLegend=legend({'x direction','y direction','z direction'},'FontSize',25,'FontWeight','bold');
hTitle =title('Force history of rear to front direction(X), left to right(Y) and from top to bottom(Z)');
% hTitle =title('Stress history of top to botom direction(Z direction)','Fontsize' ,25);
hXLabel =xlabel('number of recorded points n, loading time t=n/256 s' ,'Fontsize' ,25);
 hYLabel =ylabel('Force(N)', 'Fontsize' ,25);
 
 set(gca, 'FontName', 'Helvetica')
set([hTitle, hXLabel, hYLabel], 'FontName', 'AvantGarde')
set([hLegend, gca], 'FontSize', 20)
set([hXLabel, hYLabel], 'FontSize', 25)
set(hTitle, 'FontSize', 20, 'FontWeight' , 'bold')

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
set(gcf, 'PaperPosition', [0 0 1920 1080]); %set(gcf,'PaperPosition',[left,bottom,width,height])
