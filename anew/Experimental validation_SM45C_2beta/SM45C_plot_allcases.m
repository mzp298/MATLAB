clear;clc;close all;
load('SM45C.mat');
%% 
figure(2);
err_ben = loglog(NFben_LCF,NFben_num_LCF,'o','MarkerSize',12,'LineWidth', 3,'MarkerEdgeColor',[208 32 144]/255, 'MarkerFaceColor','none');
hold on;
err_tor = loglog(NFtor_LCF,NFtor_num_LCF,'v','MarkerSize',12,'LineWidth', 3,'MarkerEdgeColor','b', 'MarkerFaceColor','none');
err_ben = loglog(NFben_HCF,NFben_num_HCF,'o','MarkerSize',12,'LineWidth', 3,'MarkerEdgeColor',[208 32 144]/255, 'MarkerFaceColor','none');
err_tor = loglog(NFtor_HCF,NFtor_num_HCF,'v','MarkerSize',12,'LineWidth', 3,'MarkerEdgeColor','b', 'MarkerFaceColor','none');
err_ben_m = loglog(NFbenm,NF_numm,'+','MarkerSize',12,'LineWidth', 4,'MarkerEdgeColor','k', 'MarkerFaceColor','none');
err_bt_m = loglog(NF90,NF_num90,'s','MarkerSize',12,'LineWidth', 3,'MarkerEdgeColor',[139 69 19]/255, 'MarkerFaceColor','none');
err_bt_90m = loglog(NF90m,NF_num90m,'p','MarkerSize',12,'LineWidth', 3,'MarkerEdgeColor','r', 'MarkerFaceColor','none');
set(gca ,'FontSize',30);
hXLabel = xlabel('NF_{exp}','Fontsize',30, 'FontWeight' , 'bold');
hYLabel = ylabel('NF_{num}','Fontsize',30, 'FontWeight' , 'bold');
x=1e4:1000:1e8;
y0=x;
py0=loglog(x,y0,'k','LineWidth',3);
y1=2.*x;
y2=0.5.*x;
py1=loglog(x,y1, '--k','LineWidth',3);
py2=loglog(x,y2, '--k','LineWidth',3);
axis equal;
axis([1e4 1e7 1e4 1e7]);
set(gca,'xtick',[1e4 1e5 1e6 1e7]); 
set(gca,'ytick',[1e4 1e5 1e6 1e7]); 
hLegend=legend([err_ben,err_tor,err_ben_m,err_bt_m,err_bt_90m],...
    'Bending tests(R=-1)','Torsion tests(R=-1)',...
    'Bending tests(\sigma_m=196 MPa)',...
    ['Bending-torsion 90 degree ',sprintf('\n'),'out-of-phase tests'],...
    ['Bending-torsion 90 degree ',sprintf('\n'),'out-of-phase tests with mean stress'],...
   'location','bestoutside');
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
set(gcf, 'PaperPosition', [0 0 1600 800]); %set(gcf,'PaperPosition',[left,bottom,width,height])

figure(2);
saveas(gcf,'F:\Git\Anew\figures\SM45C_plot_allcases.png');

sp=actxserver('SAPI.SpVoice');
sp.Speak('done');
