clc;
clear;
close all;
format long;
% steel 10HNAP data from Jabbado thesis
%-------------------------------------------------bending--------------------------------
load('HNAP.mat');
load('gaussian.mat');
%vpa(x,289);
%vpa(weight,289);
m=1e6.*[75 ;150 ;225;300];%mean tension
m=repmat(m,1,10);
marker=['o' 'p' '^' '+'];
linestyle=['-' '--' ':' '-.' ];
color=[[208 32 144]/255; [255 140 0]/255; [255 215 0]/255; [0 139 139]/255];

NFben=[1.00E+05 2.00E+05 3.00E+05 4.00E+05 5.00E+05 6.00E+05  7.00E+05 8.00E+05 9.00E+05 1.00E+06] ;
stressben=1e6.*[311.30 289.36 276.53 267.43 260.37 254.60 249.72 245.49 241.77 238.43;...
281.77 267.18 258.64 252.58 247.89 244.05 240.80 237.99 235.51 233.29; ...
 257.98 242.47 233.40 226.96 221.97 217.89 214.44 211.46 208.82 206.47;...
251.82 224.42 208.40 197.03 188.21 181.00 174.91 169.63 164.97 160.81; ];%to get Smaxben

for  j=1:4
for  i=1:length(NFben)
    n=1;       %initial recording point
    tensor = [stressben(j,i)*sind(n*360/stepnumber)+m(j,i) 0 0 ;...
        0 0 0 ;...
        0 0 0 ];
    %---------------------to get the the first Sb-----------------------------
    run('Damiter1.m')

    while n<cycles*stepnumber
        tensor = [stressben(j,i)*sind(n*360/stepnumber)+m(j,i), 0, 0 ;...
            0, 0, 0 ;...
            0, 0, 0 ; ];
        hydro(n)=1/3*trace(tensor);
        dev1=tensor-hydro(n)*eye(3);

        dev11=dev1(1,1); dev12=dev1(1,2); dev13=dev1(1,3);%to give \dot{dev\Sigma}dt
        dev21=dev1(2,1); dev22=dev1(2,2); dev23=dev1(2,3);
        dev31=dev1(3,1); dev32=dev1(3,2); dev33=dev1(3,3);
        tensor = [stressben(j,i)*sind((n+1)*360/stepnumber)+m(j,i), 0, 0 ;...
            0, 0,  0 ;...
            0, 0,  0 ; ];
        run('Damiter2.m')
        n=n+1;
    end
    Smax_ben(j,i)=(max(Smax)-min(Smax))/2; %sqrt(2/3).*stressben
    alp_ben(j,i)=mean(alp);
	hydroplus(j,i)=mean(hydro(find(hydro>0)));
    hydrominus(j,i)=mean(hydro(find(hydro<0)));

end
parameters=[W0,b];
fun_analytical=@(parameters)parameters(1).*E.*(E+k.*nu).*parameters(2).*(parameters(2)+1).*...
    (2.*(1-alp_ben(j,:)).*(E-k).*(1+nu).*(parameters(2)-1)).^-1.*...
    (Smax_ben(j,:).^(parameters(2)+1).*(y-1/3.*lamplus.*hydroplus(j,:)).^(1-parameters(2))...
    +Smax_ben(j,:).^(parameters(2)+1).*(y-1/3.*lamminus.*hydrominus(j,:)).^(1-parameters(2))).^-1-NFben;
NFben_num(j,:)=fun_analytical(parameters)+NFben;

figure(1);
hold on;
ben_exp(j)=semilogx(NFben,Smax_ben(j,:),marker(j),'MarkerSize',12,'LineWidth', 3,'MarkerEdgeColor',color(j,:), 'MarkerFaceColor','none');
ben_num(j)=semilogx(NFben_num(j,:),Smax_ben(j,:),linestyle(j),'color',color(j,:),'LineWidth', 3);
figure(2);
err_ben_m(j) = loglog(NFben,NFben_num(j,:),marker(j),'MarkerSize',12,'LineWidth', 3,'MarkerEdgeColor',color(j,:), 'MarkerFaceColor','none');
hold on;
end

%% 

figure(1);
set(gca ,'FontSize',30);
hXLabel = xlabel('NF','Fontsize',30, 'FontWeight' , 'bold');
hYLabel = ylabel('S_{max}','Fontsize',30, 'FontWeight' , 'bold');
hLegend=legend([ben_exp(1),ben_exp(2),ben_exp(3),ben_exp(4),ben_num(1),ben_num(2),ben_num(3),ben_num(4)],...
    'Bending experimental result with \sigma_m=75 MPa','Bending experimental result with \sigma_m=150 MPa',...
    'Bending experimental result with \sigma_m=225 MPa','Bending experimental result with \sigma_m=300 MPa',...
    'Bending numerical result with \sigma_m=75 MPa','Bending numerical result with \sigma_m=150 MPa',...
    'Bending numerical result with \sigma_m=225 MPa','Bending numerical result with \sigma_m=300 MPa',...
'location','best');
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
hLegend=legend([err_ben_m(1),err_ben_m(2),err_ben_m(3),err_ben_m(4),],...
    'Bending test with mean stress \sigma_m=75 MPa on 10HNAP','Bending test with mean stress \sigma_m=150 MPa on 10HNAP',...
'Bending test with mean stress \sigma_m=225 MPa on 10HNAP','Bending test with mean stress \sigma_m=300 MPa on 10HNAP',...
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
saveas(gcf,'F:\Git\Anew\figures\10HNAP_b1D_m_sn.png');
figure(2);
saveas(gcf,'F:\Git\Anew\figures\10HNAP_b1D_m_err.png');

sp=actxserver('SAPI.SpVoice');
sp.Speak('done');


