clear;clc;
close all
format long e

fid = fopen('F:\Git\Cetim\ep_a_06\Acqui_CV.txt');
repetition=26316;
[force]=textscan(fid,'%*s%*s%s%*s%*s',repetition,'headerlines',5);
area=2.29*10*1e-6; %meter square
stress11=1000*str2double(strrep(force{1,1},',','.')).*area^-1; %Pa


ep06=plot(stress11(1:repetition),'.b','markersize',10);

grid on;
grid minor;
hTitle =title('Random loading history', 'FontName', 'AvantGarde','FontSize', 35);
axis([0 repetition -max(stress11) max(stress11)]);
hXLabel =xlabel('Number of recorded points' ,'Fontsize' ,35,'FontWeight' , 'bold');
 hYLabel =ylabel('Stress(N)', 'Fontsize' ,35,'FontWeight' , 'bold');
 
 set(gca, 'FontName', 'Helvetica', 'FontSize', 35);
set([hXLabel, hYLabel], 'FontName', 'AvantGarde');
set([hXLabel, hYLabel], 'FontSize', 35);


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
set(gcf, 'PaperPosition', [0 0 1080 600]); %set(gcf,'PaperPosition',[left,bottom,width,height])
 saveas(gcf,'F:\Git\Anew\figures\EP_a_06_random.png');
  saveas(gcf,' F:\Git\PhDreport\BEAMER\PSA 07.2017\figures\EP_a_06_random.png');
