load('HNAP.mat');%updated lambda

%------------plotting-------------------
figure(1);%----SN---
experiments_ben=semilogx(NFben,Smax_ben,'ko','MarkerSize',12,'LineWidth', 3);
hold on;
experiments_tor=semilogx(NFtor,Smax_tor,'ks','MarkerSize',12,'LineWidth', 3);
MatlabFit_ben=semilogx(NFben_num,Smax_ben,'r^','MarkerSize',12,'LineWidth', 3);
% ana_ben=plot(NFben_ana,Smax_ben,'m^','MarkerSize',12,'LineWidth', 3);
MatlabFit_tor=semilogx(NFtor_num,Smax_tor,'bv','MarkerSize',12,'LineWidth', 3);
set(gca ,'FontSize',30);
xlabel NF;
ylabel Smax;
hLegend=legend([experiments_ben,experiments_tor,MatlabFit_ben,MatlabFit_tor],...
    'Bending experiments','Torsion experiments',...
    'Bending numerical','Torsion best Fit','location','best');
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

%% 
figure(2);
err_ben = loglog(NFben,NFben_num,'o','MarkerSize',12,'LineWidth', 3,'MarkerEdgeColor','r', 'MarkerFaceColor','none');
hold on;
% err_ana_ben=loglog(NFben,NFben_ana,'m^','MarkerSize',12,'LineWidth', 3);
err_tor = loglog(NFtor,NFtor_num,'v','MarkerSize',12,'LineWidth', 3,'MarkerEdgeColor','b', 'MarkerFaceColor','none');
set(gca ,'FontSize',30);
hXLabel = xlabel('NF_{exp}','Fontsize',30, 'FontWeight' , 'bold');
hYLabel = ylabel('NF_{num}','Fontsize',30, 'FontWeight' , 'bold');
x=1e4:1000:1e8;
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
hLegend=legend([err_ben,err_tor],...
    'Bending test on 10HNAP(R=-1)','Torsion test on 10HNAP(R=-1)',...
   'location','northwest');
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
saveas(gcf,'F:\Git\Anew\figures\10HNAP_bt1D_sn.png');
figure(2);
saveas(gcf,'F:\Git\Anew\figures\10HNAP_bt1D_err.png');

sp=actxserver('SAPI.SpVoice');
sp.Speak('done');