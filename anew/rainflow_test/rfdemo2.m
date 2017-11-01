function rfdemo2(ext)
% function rfdemo2(ext)
%
% RFDEMO1 shows cycles extracted from signal
% using rainflow algoritm. Good for very long
% time signals (100 000 points).
% 
% INPUT:  ext - option; number, vector with turning
%               points or pure signal. Default ext=10000.
% 
% OUTPUT: no enable.
% 
% SYNTAX:
%         >>rfdemo2
%         >>rfdemo2(50000)
%         >>rfdemo2(my_time_signal)

% By Adam Nies³ony
% ajn@po.opole.pl

error(nargchk(0,1,nargin))

if nargin==0,
    % turning points from 10000 random numbers
    ext=sig2ext(1.7e8.*randn(1e6,1)-1e8);
elseif length(ext(:))==1,
    % turning points from n random numbers    
    ext=sig2ext(randn(1,ext));
else
    % turning points from vector ext
    ext=sig2ext(ext);
end
figure(4);
plot(ext);
grid on;
grid minor;
 hLegend=legend({'Random stress history'},'FontSize',25,'FontWeight','bold');
 set([hLegend, gca], 'FontSize', 20)
set(hLegend,'Box','on');
set(hLegend,'EdgeColor',[1 1 1]); 
hXLabel =xlabel('Number of turning points' ,'Fontsize' ,25,'FontWeight' , 'bold');
hYLabel =ylabel('Stress(Pa)', 'Fontsize' ,25,'FontWeight' , 'bold');
  set(gca, 'FontName', 'Helvetica')
set([hLegend, gca], 'FontSize', 20,'FontWeight' , 'bold');
set([hXLabel, hYLabel], 'FontSize', 20,'FontWeight' , 'bold');
set(gca, 'Box', 'off', 'TickDir', 'out', 'TickLength', [.02 .02], ...
    'XMinorTick', 'on', 'YMinorTick', 'on', 'YGrid', 'on', ...
    'XColor', [.3 .3 .3], 'YColor', [.3 .3 .3], ...
    'LineWidth', 1)
set(gcf,'color','w'); %set figure background transparent
set(gca,'color','w'); %set axis transparent
set(gcf,'outerposition',get(0,'screensize'));
set(gcf, 'PaperPositionMode', 'manual');
set(gcf, 'PaperUnits', 'points'); %[ {inches} | centimeters | normalized | points ]
set(gcf, 'PaperPosition', [0 0 1080 800]); 
 saveas(gcf,'F:\Git\Doctor_thesis_Zepeng\figures\randomdemo.png');
ext=sig2ext(ext);
rf=rainflow(ext);

figure, rfhist(rf,30,'ampl')
set(gca, 'Box', 'off', 'TickDir', 'out', 'TickLength', [.02 .02], ...
    'XMinorTick', 'on', 'YMinorTick', 'on', 'YGrid', 'on', ...
    'XColor', [.3 .3 .3], 'YColor', [.3 .3 .3], ...
    'LineWidth', 1, 'FontSize', 20,'FontWeight','bold');
set(gcf,'color','w'); %set figure background transparent
set(gca,'color','w'); %set axis transparent
set(gcf,'outerposition',get(0,'screensize'));
set(gcf, 'PaperPositionMode', 'manual');
set(gcf, 'PaperUnits', 'points'); %[ {inches} | centimeters | normalized | points ]
set(gcf, 'PaperPosition', [0 0 1080 800]); 
 saveas(gcf,'F:\Git\Doctor_thesis_Zepeng\figures\randomdemo_ampl.png');
 
figure, rfhist(rf,30,'mean')
set(gca, 'Box', 'off', 'TickDir', 'out', 'TickLength', [.02 .02], ...
    'XMinorTick', 'on', 'YMinorTick', 'on', 'YGrid', 'on', ...
    'XColor', [.3 .3 .3], 'YColor', [.3 .3 .3], ...
    'LineWidth', 1, 'FontSize', 20,'FontWeight','bold');
set(gcf,'color','w'); %set figure background transparent
set(gca,'color','w'); %set axis transparent
set(gcf,'outerposition',get(0,'screensize'));
set(gcf, 'PaperPositionMode', 'manual');
set(gcf, 'PaperUnits', 'points'); %[ {inches} | centimeters | normalized | points ]
set(gcf, 'PaperPosition', [0 0 1080 800]); 
 saveas(gcf,'F:\Git\Doctor_thesis_Zepeng\figures\randomdemo_mean.png');
 
figure, rfmatrix(rf,30,30)
set(gca, 'Box', 'off', 'TickDir', 'out', 'TickLength', [.02 .02], ...
    'XMinorTick', 'on', 'YMinorTick', 'on', 'YGrid', 'on', ...
    'XColor', [.3 .3 .3], 'YColor', [.3 .3 .3], ...
    'LineWidth', 1, 'FontSize', 20,'FontWeight','bold');
set(gcf,'color','w'); %set figure background transparent
set(gca,'color','w'); %set axis transparent
set(gcf,'outerposition',get(0,'screensize'));
set(gcf, 'PaperPositionMode', 'manual');
set(gcf, 'PaperUnits', 'points'); %[ {inches} | centimeters | normalized | points ]
set(gcf, 'PaperPosition', [0 0 1080 1080]); 
 saveas(gcf,'F:\Git\Doctor_thesis_Zepeng\figures\randomdemo_3d.png');