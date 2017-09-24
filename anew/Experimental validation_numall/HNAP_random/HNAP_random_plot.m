clear;clc;
close all
format long e
run('Parameters_HNAP')

fid = fopen('F:\Git\MATLAB\anew\Experimental validation_numall\HNAP_random\4icbmff1994Macharandomrawtimehistory.txt');
% fid = fopen('/home/ma/MATLAB/HNAP_random/4icbmff1994Macharandomrawtimehistory.txt');
%---------------------1 Numerical method with optimal time steps-----------------------------
[force]=textscan(fid,'%f',repetition,'headerlines',0);
sig_xx_original=f_sig_xx(1)*force{1,1}*1e6;%Pa
tau_xy_original=f_tau_xy(1)*force{1,1}*1e6; %Pa
% plot(force{1,1},'-k');
% grid on;
% grid minor;
% % hTitle =title('Force history on a car suspension arm in 3 directions in cartesian coordinate');
% hXLabel =xlabel('49152 random points per block with frenquency f=266.67 Hz' ,'FontWeight' , 'bold');
% hYLabel =ylabel('Signal amplitude','FontWeight' , 'bold');
% set(gca, 'FontName', 'Helvetica')
% set([hXLabel, hYLabel], 'FontName', 'AvantGarde')
% set([hXLabel, hYLabel], 'FontSize', 25)
% set(gca, 'Box', 'off', 'TickDir', 'out', 'TickLength', [.02 .02], ...
%     'XMinorTick', 'on', 'YMinorTick', 'on', 'YGrid', 'on', ...
%     'XColor', [.3 .3 .3], 'YColor', [.3 .3 .3], ...
%     'LineWidth', 1)
% set(gcf,'color','w'); %set figure background transparent
% set(gca,'color','w'); %set axis transparent
% set(gcf,'outerposition',get(0,'screensize'));
% set(gcf, 'PaperPositionMode', 'manual');
% set(gcf, 'PaperUnits', 'points'); %[ {inches} | centimeters | normalized | points ]
% set(gcf, 'PaperPosition', [0 0 1080 600]); %set(gcf,'PaperPosition',[left,bottom,width,height])
% saveas(gcf,'F:\Git\Anew\figures\10HNAPrandomblock.png');
%% 

for i=2:length(sig_xx_original)
    sig_xx(1+ari*(i-2):1+ari*(i-1))=linspace(sig_xx_original(i-1),sig_xx_original(i),ari+1);
    tau_xy(1+ari*(i-2):1+ari*(i-1))=linspace(tau_xy_original(i-1),tau_xy_original(i),ari+1);
end

% plot(sig_xx_original/1e6,'color',[238 201 0]/255)
% hold on
% plot(tau_xy_original/1e6,'color',[160 32 240]/255)
plot(sig_xx/1e6,'color',[238 201 0]/255)
hold on
plot(tau_xy/1e6,'color',[160 32 240]/255)
grid on;
grid minor;
text1=text(1.7e4,700,'$\tau_{xy}/\sigma_{xx}=0.2;\; \alpha=\pi/8$');
text2=text(0.2e4,-700,'$T\;=\;184.32\;sec$');
set(text1,'Interpreter','latex', 'Fontsize' ,35,'FontWeight' , 'bold');
set(text2,'Interpreter','latex', 'Fontsize' ,35,'FontWeight' , 'bold');
 hLegend=legend({'$\sigma_{xx}$','$\tau_{xy}$'});
 set([hLegend, gca],'FontWeight','bold')
 set(hLegend,'Interpreter','latex');
set(hLegend,'Box','on');
set(hLegend,'EdgeColor',[1 1 1]); %set the edge colour of the legend to white 
% hTitle =title('Force history on a car suspension arm in 3 directions in cartesian coordinate');
hXLabel =xlabel('49152 random points per block with frenquency f=266.67 Hz' ,'FontWeight' , 'bold');
hYLabel =ylabel('Stress(MPa)','FontWeight' , 'bold');
set(gca, 'FontName', 'Helvetica')
set([hXLabel, hYLabel], 'FontName', 'AvantGarde')
set([hLegend, gca], 'FontSize', 35)
set([hXLabel, hYLabel], 'FontSize', 25)
set(gca, 'Box', 'off', 'TickDir', 'out', 'TickLength', [.02 .02], ...
    'XMinorTick', 'on', 'YMinorTick', 'on', 'YGrid', 'on', ...
    'XColor', [.3 .3 .3], 'YColor', [.3 .3 .3], ...
    'LineWidth', 1)
set(gcf,'color','w'); %set figure background transparent
set(gca,'color','w'); %set axis transparent
set(gcf,'outerposition',get(0,'screensize'));
set(gcf, 'PaperPositionMode', 'manual');
set(gcf, 'PaperUnits', 'points'); %[ {inches} | centimeters | normalized | points ]
set(gcf, 'PaperPosition', [0 0 1080 600]); %set(gcf,'PaperPosition',[left,bottom,width,height])
% saveas(gcf,'F:\Git\Anew\figures\HNAP_random.png');
%% 
%  
% n=1;
% tensor = [sig_xx(1) tau_xy(1) 0 ;...
%     tau_xy(1) 0 0 ;...
%     0 0 0 ];
% run('Damiter1')
% while n<length(sig_xx)
%     tensor = [sig_xx(n), tau_xy(n), 0 ;...
%         tau_xy(n), 0, 0 ;...
%         0, 0, 0 ; ];
%     hydro=1/3*trace(tensor);
%     dev1=tensor-hydro*eye(3);
%     dev11=dev1(1,1); dev12=dev1(1,2); dev13=dev1(1,3);
%     dev21=dev1(2,1); dev22=dev1(2,2); dev23=dev1(2,3);
%     dev31=dev1(3,1); dev32=dev1(3,2); dev33=dev1(3,3);
%     tensor = [sig_xx(n+1), tau_xy(n+1), 0 ;...
%         tau_xy(n+1), 0,  0 ;...
%         0, 0,  0 ; ];
%     run('Damiter2')
%     n=n+1;
% end
% length(alp_ref)
% min(alp_ref)
% %% 
% 
% %----------scalar iteration
% e=1; %first D index
% j=1; %first reference alp, W, n index when iterate(after adaptation)
% D= 0;
% while D<1 %-----------the optimal time steps can be iterated with scalar
%     D= (D^(1-alp_ref(j))+(1-alp_ref(j)).*W_ref(j)/W0).^(1/(1-alp_ref(j))); % implicit
%     j=j+1;
%     if j+1>=length(alp_ref)
%         j=1;
%     end
%     e=e+1;
% end
% nF_num=(e-1)/ari;
% nF_exp=T_exp/184.32*repetition;
% xlwrite('HNAP_random.xls',nF_num,1,'G02');


