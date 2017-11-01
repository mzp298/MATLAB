clc;
clear;
close all;
format long;
% steel SM45C data from Jabbado thesis(out of phase 90)
%-------------------------------------------------bending--------------------------------
load('SM45C.mat'); %final fitting
load('gaussian.mat');

%vpa(x,289);
%vpa(weight,289);

m=0;%mean tension
NF=[29.9E3 35.7E3 50E3 73.8E3 106E3 106E3 112E3 131E3 333E3 431E3]' ;
stressben=[449E6 354E6 485E6 357E6 449E6 370E6 449E6 457E6 354E6 437E6]';%to get Smaxben
stresstor=[282E6 334E6 223E6 309E6 217E6 285E6 199E6 194E6 252E6 154E6]';%to get Smaxtor
%---------------------Numerical to get the mean value via several cycles-----------------------------


for  i=1:length(NF)
    n=1;       %initial recording point
    tensor = [stressben(i)*cosd(n*360/stepnumber)+m stresstor(i)*sind(n*360/stepnumber) 0 ;...
        stresstor(i)*sind(n*360/stepnumber) 0  0 ;...
        0 0 0 ];
		sigm=m;
    scentre=[2*sigm/3            0                0 ;...
        0             -sigm/3                0 ;...
        0                       0      -sigm/3]; %-----------ocilation center------------
    %---------------------to get the the first Sb-----------------------------
    run('Damiter1_90outofphase.m')
    while n<cycles90*stepnumber % to reach equal number of cycles as in phase
        tensor = [stressben(i)*cosd(n*360/stepnumber)+m, stresstor(i)*sind(n*360/stepnumber), 0 ;...
            stresstor(i)*sind(n*360/stepnumber), 0,  0 ;...
            0,0, 0; ];
        hydro=1/3*trace(tensor);
        dev1=tensor-hydro*eye(3)-scentre;
        dev11=dev1(1,1); dev12=dev1(1,2); dev13=dev1(1,3);%to give \dot{dev\Sigma}dt
        dev21=dev1(2,1); dev22=dev1(2,2); dev23=dev1(2,3);
        dev31=dev1(3,1); dev32=dev1(3,2); dev33=dev1(3,3);
        tensor = [stressben(i)*cosd((n+1)*360/stepnumber)+m, stresstor(i)*sind((n+1)*360/stepnumber), 0 ;...
            stresstor(i)*sind((n+1)*360/stepnumber), 0,  0 ;...
            0, 0, 0; ];
        run('Damiter2_90outofphase.m')
        n=n+1;
    end
    Smax_bt2d90(i)=max(Smax); %max sqrt of J2,a
        e=1; %first D index
    j=1; %first reference alp, W, n index when iterate(after adaptation)
    D= 1e-16;
    while D<1 %-----------the optimal time steps can be iterated with scalar
        D=D+D^alp_ref(j)*W_ref(j)/W0;
        j=j+1;
        if j>=length(alp_ref)
            j=1;
        end
        e=e+1;
    end
    NF_num(i)=e/stepnumber
end

%------------plotting-------------------
figure(1);%----SN---
experiments_bt2d90=semilogx(NF,Smax_bt2d90,'ko','MarkerSize',15,'LineWidth', 3);
hold on;
MatlabFit_bt2d90=semilogx(NF_num,Smax_bt2d90,'ms','MarkerSize',15,'LineWidth', 3);
set(gca ,'FontSize',30);
hXLabel = xlabel('NF','Fontsize',30, 'FontWeight' , 'bold');
hYLabel = ylabel('S_{a}','Fontsize',30, 'FontWeight' , 'bold');
hLegend=legend([experiments_bt2d90,MatlabFit_bt2d90],...
    ['Bending-torsion 90 degree ',sprintf('\n'),'out-of-phase experiments'],...
    ['Bending-torsion 90 degree ',sprintf('\n'),'out-of-phase numerical results'],'location','best');
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
err_bt_m = loglog(NF,NF_num,'o','MarkerSize',12,'LineWidth', 3,'MarkerEdgeColor',[139 69 19]/255, 'MarkerFaceColor','none');
set(gca ,'FontSize',30);
hXLabel = xlabel('NF_{exp}','Fontsize',30, 'FontWeight' , 'bold');
hYLabel = ylabel('NF_{num}','Fontsize',30, 'FontWeight' , 'bold');
x=1e4:1000:1e6;
y0=x;
hold on;
py0=loglog(x,y0,'k','LineWidth',3);
y1=2.*x;
y2=0.5.*x;
py1=loglog(x,y1, '--k','LineWidth',3);
py2=loglog(x,y2, '--k','LineWidth',3);
axis equal;
axis([1e4 1e6 1e4 1e6]);
set(gca,'xtick',[1e4 1e5 1e6]);
set(gca,'ytick',[1e4 1e5 1e6]);
hLegend=legend(err_bt_m,...
    ['Bending-torsion 90 degree ',sprintf('\n'),'out-of-phase test on SM45C'],'location','southeast');
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
saveas(gcf,'F:\Git\Anew\figures\SM45C_bt2D90_sn.png');
figure(2);
saveas(gcf,'F:\Git\Anew\figures\SM45C_bt2D90_err.png');


sp=actxserver('SAPI.SpVoice');
sp.Speak('done');

