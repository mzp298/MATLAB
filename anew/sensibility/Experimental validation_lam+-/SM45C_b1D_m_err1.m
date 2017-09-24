clc;
clear;
close all;
format long;
% steel SM45C data from Jabbado thesis
%-------------------------------------------------bending--------------------------------
load('SM45C.mat');
load('gaussian.mat');
%vpa(x,289);
%vpa(weight,289);

m=196e6;%mean tension
NFben=[4.50E+04
5.40E+04
7.00E+04
7.10E+04
1.10E+05
1.50E+05
2.00E+05
2.10E+05
3.00E+05
4.30E+05
6.90E+05
5.80E+05] ;
stressben=1e6.*[540
515
520
485
485
475
460
455
435
415
410
390];%to get Smaxben

for  i=1:length(NFben)
    n=1;       %initial recording point
    tensor = [stressben(i)*sind(n*360/stepnumber)+m 0 0 ;...
        0 0 0 ;...
        0 0 0 ];
    %---------------------to get the the first Sb-----------------------------
    run('Damiter1.m')
      g=1; %first reference alp, W, n index when generate  
    while n<cycles*stepnumber
        tensor = [stressben(i)*sind(n*360/stepnumber)+m, 0, 0 ;...
            0, 0, 0 ;...
            0, 0, 0 ; ];
        hydro=1/3*trace(tensor);
        dev1=tensor-hydro*eye(3);
        dev11=dev1(1,1); dev12=dev1(1,2); dev13=dev1(1,3);%to give \dot{dev\Sigma}dt
        dev21=dev1(2,1); dev22=dev1(2,2); dev23=dev1(2,3);
        dev31=dev1(3,1); dev32=dev1(3,2); dev33=dev1(3,3);
           tensor = [stressben(i)*sind((n+1)*360/stepnumber)+m, 0, 0 ;...
            0, 0,  0 ;...
            0, 0,  0 ; ];
       run('Damiter2.m')
                n=n+1;
    end
    Smax_ben_m(i)=(max(Smax)-min(Smax))/2; %sqrt(2/3).*stressben
    alp_ben_m(i)=mean(alp);
    hydroplus(i)=mean(hydro(find(hydro>0)));
    hydrominus(i)=mean(hydro(find(hydro<0)));
end
parameters=[W0,b];
fun_analytical=@(parameters)parameters(1).*E.*(E+k.*nu).*parameters(2).*(parameters(2)+1).*...
    (2.*(1-alp_ben_m).*(E-k).*(1+nu).*(parameters(2)-1)).^-1.*...
    (Smax_ben_m.^(parameters(2)+1).*(y-1/3.*lamplus.*hydroplus).^(1-parameters(2))...
    +Smax_ben_m.^(parameters(2)+1).*(y-1/3.*lamminus.*hydrominus).^(1-parameters(2))).^-1-NF';
NF_num=fun_analytical(parameters)+NF'


%% 
%--------use average smin to get NF--------
figure(1);
smax_exp=semilogx(NFben(1:i),Smax_ben_m(1:i),'o','MarkerSize',12,'LineWidth', 3,'MarkerEdgeColor',[208 32 144]/255, 'MarkerFaceColor','none');
hold on;
smax_num=semilogx(NF_num(1:i),Smax_ben_m(1:i),'-r','LineWidth', 3);

set(gca ,'FontSize',30);
hXLabel = xlabel('NF_{num}','Fontsize',30, 'FontWeight' , 'bold');
hYLabel = ylabel('S_{max}','Fontsize',30, 'FontWeight' , 'bold');
hLegend=legend([smax_exp,smax_num],...
    'Bending experimental result with \sigma_m=196 MPa','Bending numerical result with \sigma_m=196 MPa','location','best');
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
err_ben_m = loglog(NFben,NF_num,'o','MarkerSize',12,'LineWidth', 3,'MarkerEdgeColor',[208 32 144]/255, 'MarkerFaceColor','none');
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
    'Bending test on SM45(\sigma_m=196 MPa)','location','northwest');
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
saveas(gcf,'F:\Git\Anew\figures\SM45C_b1D_m_sn.png');
figure(2);
saveas(gcf,'F:\Git\Anew\figures\SM45C_b1D_m_err.png');

sp=actxserver('SAPI.SpVoice');
sp.Speak('done');
