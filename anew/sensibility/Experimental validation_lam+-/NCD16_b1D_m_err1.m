clc;
clear;
close all;
format long;
% steel 30 NCD 16 data from Jabbado thesis
%-------------------------------------------------bending--------------------------------
load('NCD16.mat','lam','W0','b','a','E','k','nu','y','stepnumber','cycles','delta_alp'); %final fitting
load('gaussian.mat');
%vpa(x,289);
%vpa(weight,289);

m=1e6.*[450 450 290 290];%mean tension
NFben=[51000 140000 120000 250000] ;
stressben=1e6.*[640 620 695 660];%to get Smaxben

for  i=1:length(NFben)
    n=1;       %initial recording point
    tensor = [stressben(i)*sind(n*360/stepnumber)+m(i) 0 0 ;...
        0 0 0 ;...
        0 0 0 ];
    run('Damiter1.m')
    g=1; %first reference alp, W, n index when generate
    while n<cycles*stepnumber
        tensor = [stressben(i)*sind(n*360/stepnumber)+m(i), 0, 0 ;...
            0, 0, 0 ;...
            0, 0, 0 ; ];
        hydro(n)=1/3*trace(tensor);
        dev1=tensor-hydro(n)*eye(3);

        dev11=dev1(1,1); dev12=dev1(1,2); dev13=dev1(1,3);%to give \dot{dev\Sigma}dt
        dev21=dev1(2,1); dev22=dev1(2,2); dev23=dev1(2,3);
        dev31=dev1(3,1); dev32=dev1(3,2); dev33=dev1(3,3);
        tensor = [stressben(i)*sind((n+1)*360/stepnumber)+m(i), 0, 0 ;...
            0, 0,  0 ;...
            0, 0,  0 ; ];
        run('Damiter2.m')
        n=n+1;
    end
    Smax_ben(i)=(max(Smax)-min(Smax))/2; %sqrt(2/3).*stressben
    %     Wcyc1(i)=4*(E-k)*(1+nu)*(b-1)/(E*(E+k*nu)*b*(b+1)).*Smax_ben(i).^(b+1).*(y-1/3*lam*m(i)).^(1-b) ;
    alp_ben(i)=mean(alp);
    hydroplus(i)=mean(hydro(find(hydro>0)));
    hydrominus(i)=mean(hydro(find(hydro<0)));
end
parameters=[W0,b];
fun_analytical=@(parameters)parameters(1).*E.*(E+k.*nu).*parameters(2).*(parameters(2)+1).*...
    (2.*(1-alp_ben).*(E-k).*(1+nu).*(parameters(2)-1)).^-1.*...
    (Smax_ben.^(parameters(2)+1).*(y-1/3.*lamplus.*hydroplus).^(1-parameters(2))...
    +Smax_ben.^(parameters(2)+1).*(y-1/3.*lamminus.*hydrominus).^(1-parameters(2))).^-1-NF';
NF_num=fun_analytical(parameters)+NF'

%------------plotting-------------------
figure(1);%----SN---
experiments_ben=plot(NFben,Smax_ben,'ko','MarkerSize',12,'LineWidth', 3);
hold on;
MatlabFit_ben=plot(NFben_num,Smax_ben,'r^','MarkerSize',12,'LineWidth', 3);
xlabel NF;
ylabel Smax;
hLegend=legend([experiments_ben,MatlabFit_ben],...
    'Bending-torsion experiments',...
    'Bending-torsion best Fit','location','best');
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
err_ben_m = loglog(NFben,NFben_num,'o','MarkerSize',12,'LineWidth', 3,'MarkerEdgeColor',[208 32 144]/255, 'MarkerFaceColor','none');
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
hLegend=legend([err_ben_m],...
    'Bending test on 30NCD16 steel with mean stress','location','northwest');
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
saveas(gcf,'F:\Git\Anew\figures\NCD16_b1D_m_sn.png');
figure(2);
saveas(gcf,'F:\Git\Anew\figures\NCD16_b1D_m_err.png');

sp=actxserver('SAPI.SpVoice');
sp.Speak('done');
