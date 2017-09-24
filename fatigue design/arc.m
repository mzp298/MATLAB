clear;clc;
syms tau sig;
sref=397;
f=460;
tref=258;
t=292;

a=(tref-sref/sqrt(3))/(sref/3);
b=tref;
ap=3*(tref/f-1/2);
bp=tref;
ad=3*(tref/sref-1/2);
bd=tref;
y=20;
lg=2.5;
lgp=2.5;
lgd=2.5;

cross=(tau/tref)^2+(2*sref/(sqrt(3)*tref)-1)*(sig/sref)^2+(2-2*sref/(sqrt(3)*tref))*(sig/sref)-1;%classiclal crossland
papa= sqrt(sig^2/3+tau^2)+ap*sig/3-bp;%classiclal papa
dangvan=sqrt(sig^2/4+tau^2)+ad*sig/3-bd;%classiclal dangvan

papaf=(tau/tref)^2+(2*f/(sqrt(3)*tref)-1)*(sig/f)^2+(2-2*f/(sqrt(3)*tref))*(sig/f)-1;%papadopoulos criterion using f_{-1} instead of s_{-1}

crossgrad=sqrt(sig^2/3+tau^2)+a*sig/3-lg*(sqrt(sig^2/3+tau^2)/y+a*sig/(3*y))-b;
papagrad= sqrt(sig^2/3+tau^2)+ap*sig/3-lgp*(sqrt(sig^2/3+tau^2)/y+ap*sig/(3*y))-bp;
dangvangrad=sqrt(sig^2/4+tau^2)+ad*sig/3-lgd*(sqrt(sig^2/4+tau^2)/y+ad*sig/(3*y))-bd;

hold on	
Crossland_Classical = set(ezplot(cross,[0,500,0,350]),'Color',[193 205 193]/255,'LineStyle', '-', 'LineWidth', 4);
Papadopoulos_tf = set(ezplot(papaf,[0,500,0,350]),'Color','k','LineStyle', '--', 'LineWidth',4);
Crossland_Gradient = set(ezplot(crossgrad,[0,500,0,350]),'Color',[238 99 99]/255,'LineStyle', '-', 'LineWidth', 4);
Papadopoulos_Gradient = set(ezplot(papagrad,[0,500,0,350]) ,'Color',[160 32 240]/255,'LineStyle', ' :', 'LineWidth', 4);
DangVan_Gradient = set(ezplot(dangvangrad,[0,500,0,350]),'Color',[238 180 34]/255,'LineStyle', '-.', 'LineWidth', 4);
sigma = [0 155 260 360 f sref];
tau     = [t 270 235 175 0 0];
exp = plot(sigma, tau, 'o','MarkerSize',20, 'MarkerFaceColor','k', 'MarkerEdgeColor',[1 0.5 0]);
hold off;

grid on;
set(gca ,'FontSize',30);
hTitle = title({['SAE 4340 steel'],[ 'Fully reversed bending and torsion'],[ '(data from Findley)']; });
set(hTitle, 'FontSize', 23, 'FontWeight' , 'bold');
hXLabel =  xlabel('\sigma_a(MPa)','Fontsize',30);
hYLabel =  ylabel('\tau_a(Mpa)','Fontsize',30);
hLegend = legend('Crossland\_Classical','Papadopoulos\_based on (t,f)','Crossland\_Gradient',...
    'Papadopoulos\_Gradient','DangVan\_Gradient','Experiments','Location','bestoutside');
set([hLegend, gca], 'FontSize', 20)
set(hLegend,'Box','off');
text(340,22,'S_{ref}','Fontsize',30);
text(423,22,'f_{-1}','Fontsize',30);
text(8,t+17,'t','Fontsize',30);
text(8,tref-22,'t_{ref}','Fontsize',30);

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
set(gcf, 'PaperPosition', [0 0 1080 565]); %set(gcf,'PaperPosition',[left,bottom,width,height])
 saveas(gcf,'4340.png');
  saveas(gcf,'4340.fig');


