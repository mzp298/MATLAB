clear;clc;close all;
Smax=1:1e4:1.71e8;
b=1.5;  
y=230e6;
yield=y;
a=0.1;
sequence=((Smax.*yield.^-1).*(1-Smax.*yield.^-1).^-1).^b;
sequence2=((Smax.*yield.^-1).*(1-Smax.*yield.^-1).^-1);
alp=1-a*sequence;
hold on
plot(Smax,a*sequence,'r','LineWidth', 5);
plot(Smax,a*sequence2,'-.b','LineWidth', 5);
grid on;
hXLabel = xlabel('S_{max}');
hYLabel =ylabel('(1-\alpha)');

hLegend=legend( '(1-\alpha) with magnification power \beta=1.1',...
    '(1-\alpha) without magnification power','Location','Best');
set([hLegend, gca], 'FontSize', 35)
set(hLegend,'Box','on');
set(hLegend,'EdgeColor',[1 1 1]); %set the edge colour of the legend to white 
% Adjust font
set(gca, 'FontName', 'Helvetica')
set([ hXLabel, hYLabel], 'FontName', 'AvantGarde')
set([hXLabel, hYLabel], 'FontSize',35)
% Adjust axes properties
set(gcf,'color','w'); %set figure background transparent
set(gca,'color','w'); %set axis transparent
% Maximize print figure
set(gcf,'outerposition',get(0,'screensize'));
set(gcf, 'PaperPositionMode', 'manual');
set(gcf, 'PaperUnits', 'points'); %[ {inches} | centimeters | normalized | points ]
set(gcf, 'PaperPosition', [0 0 1200 700]); %set(gcf,'PaperPosition',[left,bottom,width,height])
 saveas(gcf,'F:\Git\Anew\figures\alp_Smax_fb.png');
 saveas(gcf,'F:\Git\PhDreport\BEAMER\UTC2017\figures\alp_Smax_fb.png');