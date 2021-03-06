clc;
clear;
close all;
format long;
% steel 10HNAP data from Jabbado thesis
%-------------------------------------------------bending--------------------------------
load('HNAP.mat');
%  cycles=3; %--------to reach accomondation----
load('gaussian.mat');
%vpa(x,289);
%vpa(weight,289);
m=1e6.*[75 ;150 ;225;300];%mean tension
m=repmat(m,1,10);
marker=['o' 'p' 's' '^'];
linestyle=['-' '--' ':' '-.'];
color=[[208 32 144]/255; [255 140 0]/255; [255 215 0]/255; [0 139 139]/255];

NFbenm=[1.00E+05 2.00E+05 3.00E+05 4.00E+05 5.00E+05 6.00E+05  7.00E+05 8.00E+05 9.00E+05 1.00E+06] ;
stressbenm=1e6.*[311.30 289.36 276.53 267.43 260.37 254.60 249.72 245.49 241.77 238.43;...
    281.77 267.18 258.64 252.58 247.89 244.05 240.80 237.99 235.51 233.29; ...
    257.98 242.47 233.40 226.96 221.97 217.89 214.44 211.46 208.82 206.47;...
    251.82 224.42 208.40 197.03 188.21 181.00 174.91 169.63 164.97 160.81; ];%to get Smaxben

for  j=1:4 %-----------4 mean stress-----------------
    sigm=m(j);
    for  i=1:length(NFbenm)
        n=1;       %initial recording point
        tensor = [stressbenm(j,i)*sind(n*360/stepnumber)+m(j,i) 0 0 ;...
            0 0 0 ;...
            0 0 0 ];
        run('Damiter1.m')
        while n<cycles*stepnumber
            tensor = [stressbenm(j,i)*sind(n*360/stepnumber)+m(j,i), 0, 0 ;...
                0, 0, 0 ;...
                0, 0, 0 ; ];
            hydro(n)=1/3*trace(tensor);
            dev1=tensor-hydro(n)*eye(3)-scentre;
            dev11=dev1(1,1); dev12=dev1(1,2); dev13=dev1(1,3);
            dev21=dev1(2,1); dev22=dev1(2,2); dev23=dev1(2,3);
            dev31=dev1(3,1); dev32=dev1(3,2); dev33=dev1(3,3);
            tensor = [stressbenm(j,i)*sind((n+1)*360/stepnumber)+m(j,i), 0, 0 ;...
                0, 0,  0 ;...
                0, 0,  0 ; ];
            run('Damiter2.m')
            n=n+1;
        end
        Smax_ben(j,i)=max(Smax); %max sqrt of J2,a
        
        alp_ben(j,i)=mean(alp_ref);
        hydroplus(j,i)=1/sqrt(2)*(max(hydro)-mean(hydro))+mean(hydro);
        hydrominus(j,i)=1/sqrt(2)*(min(hydro)-mean(hydro))+mean(hydro);
        
        e=1; %first D index
        r=1; %first reference alp, W, n index when iterate(after adaptation)
        D= 1e-16;
        while D<1 %-----------the optimal time steps can be iterated with scalar
            D=D+D^alp_ref(r)*W_ref(r)/W0;
            r=r+1;
            if r+1>=length(alp_ref)
                r=1;
            end
            e=e+1;
        end
        NFben_m_num(j,i)=e/stepnumber
    end
end

save('HNAP.mat','NFbenm','NFben_m_num','-append');
% % -------------fit directly on 1d with mean stress(1 lam)
NFben_tensor=repmat(NFbenm,4,1);
%% 

for  j=1:4
    figure(1);
    ben_exp(j)=semilogx(NFbenm,Smax_ben(j,:),marker(j),'MarkerSize',12,'LineWidth', 3,'MarkerEdgeColor',color(j,:), 'MarkerFaceColor','none');
    hold on;
    ben_num(j)=semilogx(NFben_m_num(j,:),Smax_ben(j,:),linestyle(j),'color',color(j,:),'LineWidth', 3);
    figure(2);
    err_ben_m(j) = loglog(NFbenm,NFben_m_num(j,:),marker(j),'MarkerSize',12,'LineWidth', 3,'MarkerEdgeColor',color(j,:), 'MarkerFaceColor','none');
    hold on;
    figure(3);
    hold on;
    Sa(j) = plot(Smax_ben(j,:),marker(j),'MarkerSize',12,'LineWidth', 3,'MarkerEdgeColor',color(j,:), 'MarkerFaceColor','none');
    figure(4);
    hold on;
    hydroplusplot(j) = plot(hydroplus(j,:),marker(j),'MarkerSize',12,'LineWidth', 3,'MarkerEdgeColor',color(j,:), 'MarkerFaceColor',color(j,:));
    hydrominusplot(j) = plot(hydrominus(j,:),marker(j),'MarkerSize',12,'LineWidth', 3,'MarkerEdgeColor',color(j,:), 'MarkerFaceColor','none');
end


%%

figure(3);
grid on;
set(gca ,'FontSize',30);
hXLabel = xlabel('Test number','Fontsize',30, 'FontWeight' , 'bold');
hYLabel = ylabel('S_{a}','Fontsize',30, 'FontWeight' , 'bold');
hLegend=legend([Sa(1),Sa(2),Sa(3),Sa(4)],...
    'S_a with \sigma_m=75 MPa','S_a with \sigma_m=150 MPa',...
    'S_a with \sigma_m=225 MPa','S_a with \sigma_m=300 MPa',...
    'location','northeast');
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
set(gcf,'outerposition',get(0,'screensize'));
set(gcf, 'PaperPositionMode', 'manual');
set(gcf, 'PaperUnits', 'points'); %[ {inches} | centimeters | normalized | points ]
set(gcf, 'PaperPosition', [0 0 800 800]); %set(gcf,'PaperPosition',[left,bottom,width,height])

figure(4);
grid on;
set(gca ,'FontSize',30);
hXLabel = xlabel('Test number','Fontsize',30, 'FontWeight' , 'bold');
hYLabel = ylabel('Hydro+(solid) and hydro-(void)','Fontsize',30, 'FontWeight' , 'bold');
hLegend=legend([hydroplusplot(1),hydroplusplot(2),hydroplusplot(3),hydroplusplot(4),...
    hydrominusplot(1),hydrominusplot(2),hydrominusplot(3),hydrominusplot(4)],...
    'hydro+ with \sigma_m=75 MPa','hydro+ with \sigma_m=150 MPa',...
    'hydro+ with \sigma_m=225 MPa','hydro+ with \sigma_m=300 MPa',...
    'hydro- with \sigma_m=75 MPa','hydro- with \sigma_m=150 MPa',...
    'hydro- with \sigma_m=225 MPa','hydro- with \sigma_m=300 MPa',...
    'location','northeast');
set(hLegend, 'FontSize', 18);
set(hLegend,'Box','on');
set(hLegend,'EdgeColor',[1 1 1]); %set the edge colour of the legend to white
set(gcf,'color','w'); %set figure background transparent
set(gca,'color','w'); %set axis transparent
set(gcf,'outerposition',get(0,'screensize'));
set(gcf, 'PaperPositionMode', 'manual');
set(gcf, 'PaperUnits', 'points'); %[ {inches} | centimeters | normalized | points ]
set(gcf, 'PaperPosition', [0 0 800 800]); %set(gcf,'PaperPosition',[left,bottom,width,height])


figure(1);
set(gca ,'FontSize',30);
hXLabel = xlabel('NF','Fontsize',30, 'FontWeight' , 'bold');
hYLabel = ylabel('S_{max}','Fontsize',30, 'FontWeight' , 'bold');
hLegend=legend([ben_exp(1),ben_exp(2),ben_exp(3),ben_exp(4),ben_num(1),ben_num(2),ben_num(3),ben_num(4)],...
    'Bending_{exp} with \sigma_m=75 MPa','Bending_{exp} with \sigma_m=150 MPa',...
    'Bending_{exp} with \sigma_m=225 MPa','Bending_{exp} with \sigma_m=300 MPa',...
    'Bending_{fitting} with \sigma_m=75 MPa','Bending_{fitting} with \sigma_m=150 MPa',...
    'Bending_{fitting} with \sigma_m=225 MPa','Bending_{fitting} with \sigma_m=300 MPa',...
    'location','best');
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
set(gcf, 'PaperPosition', [0 0 800 800]); %set(gcf,'PaperPosition',[left,bottom,width,height])

figure(2);
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
hLegend=legend([err_ben_m(1),err_ben_m(2),err_ben_m(3),err_ben_m(4)],...
    'Bending with \sigma_m=75 MPa','Bending with \sigma_m=150 MPa',...
    'Bending with \sigma_m=225 MPa','Bending with \sigma_m=300 MPa',...
    'location','northwest');
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
set(gcf, 'PaperPosition', [0 0 800 800]); %set(gcf,'PaperPosition',[left,bottom,width,height])

figure(1);
saveas(gcf,'F:\Git\Anew\figures\10HNAP_b1D_m_sn.png');
figure(2);
saveas(gcf,'F:\Git\Anew\figures\10HNAP_b1D_m_err.png');
figure(3);
saveas(gcf,'F:\Git\Anew\figures\10HNAP_b1D_m_Smax.png');
figure(4);
saveas(gcf,'F:\Git\Anew\figures\10HNAP_b1D_m_hydro.png');

sp=actxserver('SAPI.SpVoice');
sp.Speak('done');


