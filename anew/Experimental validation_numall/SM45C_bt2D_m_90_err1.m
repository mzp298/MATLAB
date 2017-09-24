clc;
clear;
close all;
format long;
% steel SM45C data from Jabbado thesis(out of phase 90 degrees)
%-------------------------------------------------bending--------------------------------
load('SM45C.mat'); %final fitting
load('gaussian.mat');
%vpa(x,289);
%vpa(weight,289);

m=196e6;%mean tension
NF=[53E3 59.2E3 70.1E3 86.3E3 89.9E3 92.1E3 102E3 135E3 351E3 394E3]' ;
stressben=[441E6 286E6 464E6 473E6 173E6 403E6 437E6 167E6 357E6 182E6]';%to get Smaxben
stresstor=[215E6 309E6 155E6 136E6 334E6 209E6 177E6 321E6 179E6 274E6]';%to get Smaxtor
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
    run('Damiter1.m')
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
        run('Damiter2.m')
        n=n+1;
    end
    Smax_bt2dm90(i)=max(Smax); %max sqrt of J2,a
    alp_bt2dm90(i)=mean(alp_ref);
    hydroplus(i)=1/sqrt(2)*(max(hydro)-mean(hydro))+mean(hydro);
    hydrominus(i)=1/sqrt(2)*(min(hydro)-mean(hydro))+mean(hydro);
end
%Smax_bt2dm90=sqrt(stressben.^2+stresstor.^2)'; %out of phase Smax
fun_analytical=@(parameters)parameters(1).*E.*(E+k.*nu).*parameters(2).*(parameters(2)+1).*...
    (2.*(1-alp_bt2dm90).*(E-k).*(1+nu).*(parameters(2)-1)).^-1.*...
    (Smax_bt2dm90.^(parameters(2)+1).*(y-lamplus.*hydroplus).^(1-parameters(2))...
    +Smax_bt2dm90.^(parameters(2)+1).*(y-lamminus.*hydrominus).^(1-parameters(2))).^-1-NF';
parameters=[W0,b];
NF_num=fun_analytical(parameters)+NF';


%------------plotting-------------------
figure(1);%----SN---
experiments_bt2dm90=plot(NF,Smax_bt2dm90,'ko','MarkerSize',12,'LineWidth', 3);
hold on;
MatlabFit_bt2dm90=plot(NF_num,Smax_bt2dm90,'--m','MarkerSize',12,'LineWidth', 3);
xlabel NF;
ylabel Smax;
hLegend=legend([experiments_bt2dm90,MatlabFit_bt2dm90],...
    ['Bending-torsion 90 degree ',sprintf('\n'),'out-of-phase test on SM45C',sprintf('\n'),'with mean stress'],...
    ['Bending-torsion 90 degree ',sprintf('\n'),'out-of-phase test on SM45C',sprintf('\n'),'with mean stress fitting'],'location','best');
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
    ['Bending-torsion 90 degree ',sprintf('\n'),'out-of-phase test on SM45C',sprintf('\n'),'with mean stress'],'location','northwest');
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
saveas(gcf,'F:\Git\Anew\figures\SM45C_bt2D90_m_sn.png');
figure(2);
saveas(gcf,'F:\Git\Anew\figures\SM45C_bt2D90_m_err.png');

sp=actxserver('SAPI.SpVoice');
sp.Speak('done');

