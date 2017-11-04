clear;clc;close all;
load('HNAP.mat');
marker=['o' 'p' 's' '^'];
color=[[208 32 144]/255; [255 140 0]/255; [255 215 0]/255; [0 139 139]/255];
figure(2);%------error bar---------
for  j=1:4
 figure(2);
    err_ben_m(j) = loglog(NFbenm,NFben_m_num(j,:),marker(j),'MarkerSize',12,'LineWidth', 3,'MarkerEdgeColor',color(j,:), 'MarkerFaceColor','none');
    hold on;
end
err_ben = loglog(NFben,NFben_num,'o','MarkerSize',12,'LineWidth', 3,'MarkerEdgeColor','r', 'MarkerFaceColor','none');
err_tor = loglog(NFtor,NFtor_num,'v','MarkerSize',12,'LineWidth', 3,'MarkerEdgeColor','b', 'MarkerFaceColor','none');
set(gca ,'FontSize',30);
hXLabel = xlabel('NF_{exp}','Fontsize',30, 'FontWeight' , 'bold');
hYLabel = ylabel('NF_{num}','Fontsize',30, 'FontWeight' , 'bold');
x=1e4:1e4:1e8;
y0=x;
hold on;
py0=loglog(x,y0,'k','LineWidth',3);
y1=2.*x;
y2=0.5.*x;
py1=loglog(x,y1, '--k','LineWidth',3);
py2=loglog(x,y2, '--k','LineWidth',3);
axis equal;
axis([1e4 1e7 1e4 1e7]);
set(gca,'xtick',[1e4 1e5 1e6 1e7]);
set(gca,'ytick',[1e4 1e5 1e6 1e7]);
hLegend=legend([err_ben,err_tor,err_ben_m(1),err_ben_m(2),err_ben_m(3),err_ben_m(4),],...
    'Bending test on 10HNAP(R=-1)','Torsion test on 10HNAP(R=-1)',...
    'Bending with \sigma_m=75 MPa','Bending with \sigma_m=150 MPa',...
    'Bending with \sigma_m=225 MPa','Bending with \sigma_m=300 MPa',...
    'location','bestoutside');
set(hLegend, 'FontSize', 18);
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
set(gcf, 'PaperPosition', [0 0 1600 800]); %set(gcf,'PaperPosition',[left,bottom,width,height])

figure(2);
saveas(gcf,'F:\Git\Anew\figures\10HNAP_plot_allcases.png');


sp=actxserver('SAPI.SpVoice');
sp.Speak('done');