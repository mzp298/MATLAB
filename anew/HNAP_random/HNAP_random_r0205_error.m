clear;clc;close all;
Texp02= xlsread('1HNAP_random.xls',1,'E2:E22')
Texp05= xlsread('1HNAP_random.xls',1,'E23:E36')
Tnum02= xlsread('1HNAP_random.xls',1,'G2:G22')
Tnum05= xlsread('1HNAP_random.xls',1,'G23:G36')
JBexp02= xlsread('1HNAP_random.xls',1,'E39:E59')
JBexp05= xlsread('1HNAP_random.xls',1,'E60:E73')
JBnum02= xlsread('1HNAP_random.xls',1,'G39:G59')
JBnum05= xlsread('1HNAP_random.xls',1,'G60:G73')
%------------plotting-------------------
figure(1);%------error bar of r=0.2---------
err02 = loglog(Texp02,Tnum02,'o','MarkerSize',12,'LineWidth', 3,'MarkerEdgeColor',[208 32 144]/255, 'MarkerFaceColor','none');
hold on;
JBerr02 = loglog(JBexp02,JBnum02,'s','MarkerSize',12,'LineWidth', 3,'MarkerEdgeColor','k', 'MarkerFaceColor','none');
set(gca ,'FontSize',30);
hXLabel = xlabel('T_{exp}','Fontsize',30, 'FontWeight' , 'bold');
hYLabel = ylabel('T_{num}','Fontsize',30, 'FontWeight' , 'bold');
x=1e3:1000:1e7;
y0=x;
py0=loglog(x,y0,'k','LineWidth',3);
y1=2.*x;
y2=0.5.*x;
py1=loglog(x,y1, '--k','LineWidth',3);
py2=loglog(x,y2, '--k','LineWidth',3);
axis equal;
axis([1e3 1e7 1e3 1e7]);
set(gca,'xtick',[1e3 1e4 1e5 1e6 1e7]);
set(gca,'ytick',[1e3 1e4 1e5 1e6 1e7]);
hLegend=legend([err02,JBerr02],...
    ['Proposed model ($\alpha=\pi/8, r=0.2$)'],...
    ['Jabbado model ($\alpha=\pi/8, r=0.2$)'],...
    'location','northwest');
set(hLegend,'Interpreter','latex');
set(hLegend, 'FontSize', 28);
set(hLegend,'Box','on');
set(hLegend,'EdgeColor',[1 1 1]); %set the edge colour of the legend to white
% Adjust font
set(gca, 'FontName', 'Helvetica')
% Adjust axes properties
set(gca, 'Box', 'on', 'TickDir', 'out', 'TickLength', [.02 .02], ...
    'XMinorTick', 'on', 'YMinorTick', 'on', 'YGrid', 'on', 'XGrid', 'on',...
    'XColor', [.3 .3 .3], 'YColor', [.3 .3 .3], ...
    'LineWidth', 1)
set(gcf,'color','w'); %set figure background transparent
set(gca,'color','w'); %set axis transparent
% Maximize print figure
set(gcf,'outerposition',get(0,'screensize'));
set(gcf, 'PaperPositionMode', 'manual');
set(gcf, 'PaperUnits', 'points'); %[ {inches} | centimeters | normalized | points ]
set(gcf, 'PaperPosition', [0 0 800 800]); %set(gcf,'PaperPosition',[left,bottom,width,height])

%------------plotting-------------------
figure(2);%------error bar of r=0.5---------
err05 = loglog(Texp05,Tnum05,'o','MarkerSize',12,'LineWidth', 3,'MarkerEdgeColor',[208 32 144]/255, 'MarkerFaceColor','none');
hold on;
JBerr05 = loglog(JBexp05,JBnum05,'s','MarkerSize',12,'LineWidth', 3,'MarkerEdgeColor','k', 'MarkerFaceColor','none');
set(gca ,'FontSize',30);
hXLabel = xlabel('T_{exp}','Fontsize',30, 'FontWeight' , 'bold');
hYLabel = ylabel('T_{num}','Fontsize',30, 'FontWeight' , 'bold');
x=1e3:1000:1e7;
y0=x;
py0=loglog(x,y0,'k','LineWidth',3);
y1=2.*x;
y2=0.5.*x;
py1=loglog(x,y1, '--k','LineWidth',3);
py2=loglog(x,y2, '--k','LineWidth',3);
axis equal;
axis([1e3 1e7 1e3 1e7]);
set(gca,'xtick',[1e3 1e4 1e5 1e6 1e7]);
set(gca,'ytick',[1e3 1e4 1e5 1e6 1e7]);
hLegend=legend([err05,JBerr05],...
    ['Proposed model ($\alpha=\pi/4, r=0.5$)'],...
    ['Jabbado model ($\alpha=\pi/4, r=0.5$)'],...
    'location','northwest');
set(hLegend,'Interpreter','latex');
set(hLegend, 'FontSize', 28);
set(hLegend,'Box','on');
set(hLegend,'EdgeColor',[1 1 1]); %set the edge colour of the legend to white
% Adjust font
set(gca, 'FontName', 'Helvetica')
% Adjust axes properties
set(gca, 'Box', 'on', 'TickDir', 'out', 'TickLength', [.02 .02], ...
    'XMinorTick', 'on', 'YMinorTick', 'on', 'YGrid', 'on', 'XGrid', 'on',...
    'XColor', [.3 .3 .3], 'YColor', [.3 .3 .3], ...
    'LineWidth', 1)
set(gcf,'color','w'); %set figure background transparent
set(gca,'color','w'); %set axis transparent
% Maximize print figure
set(gcf,'outerposition',get(0,'screensize'));
set(gcf, 'PaperPositionMode', 'manual');
set(gcf, 'PaperUnits', 'points'); %[ {inches} | centimeters | normalized | points ]
set(gcf, 'PaperPosition', [0 0 800 800]); %set(gcf,'PaperPosition',[left,bottom,width,height])


figure(1);
saveas(gcf,'F:\Git\Doctor_thesis_Zepeng\figures\HNAP_random_r02_error.png');
figure(2);
saveas(gcf,'F:\Git\Doctor_thesis_Zepeng\figures\HNAP_random_r05_error.png');

sp=actxserver('SAPI.SpVoice');
sp.Speak('done');