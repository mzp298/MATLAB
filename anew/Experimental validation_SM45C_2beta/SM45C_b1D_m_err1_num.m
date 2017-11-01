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
NFben=[43043.40644
55522.84924
74724.73326
75361.55953
110407.2374
146089.5545
194951.3693
212217.5648
297989.8518
440285.7643
678727.3216
597603.0102
] ;
stressben=1e6.*[541.48629
511.47186
514.35786
493.57864
490.69264
471.06782
455.48341
452.0202
430.08658
409.30736
407.57576
386.79654
];%to get Smaxben


for  i=1:length(NFben)
    n=1;
    tensor = [stressben(i)*sind(n*360/stepnumber)+m 0 0 ;...
        0 0 0 ;...
        0 0 0 ];
    sigm=m;
    scentre=[2*sigm/3            0                0 ;...
        0             -sigm/3                0 ;...
        0                       0      -sigm/3];
    run('Damiter1.m')
    while n<cycles*stepnumber
        tensor = [stressben(i)*sind(n*360/stepnumber)+m, 0, 0 ;...
            0, 0, 0 ;...
            0, 0, 0 ; ];
        hydro=1/3*trace(tensor);
        dev1=tensor-hydro*eye(3)-scentre;
        dev11=dev1(1,1); dev12=dev1(1,2); dev13=dev1(1,3);
        dev21=dev1(2,1); dev22=dev1(2,2); dev23=dev1(2,3);
        dev31=dev1(3,1); dev32=dev1(3,2); dev33=dev1(3,3);
        tensor = [stressben(i)*sind((n+1)*360/stepnumber)+m, 0, 0 ;...
            0, 0,  0 ;...
            0, 0,  0 ; ];
        run('Damiter2.m')
        n=n+1;
    end
    Smax_ben_m(i)=max(Smax);
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

figure(1);
smax_exp=semilogx(NFben,Smax_ben_m,'ko','MarkerSize',15,'LineWidth', 3,'MarkerFaceColor','none');
hold on;
smax_num=semilogx(NF_num,Smax_ben_m,'b^','MarkerSize',15,'LineWidth', 3);
set(gca ,'FontSize',30);
hXLabel = xlabel('NF','Fontsize',30, 'FontWeight' , 'bold');
hYLabel = ylabel('S_{a}','Fontsize',30, 'FontWeight' , 'bold');
hLegend=legend([smax_exp,smax_num],...
    ['Bending experimental result ',sprintf('\n'),'with \sigma_m=196 MPa'],...
    ['Bending numerical result ',sprintf('\n'),'with \sigma_m=196 MPa'],'location','best');
set(hLegend, 'FontSize', 28);
set(hLegend,'Box','on');
set(hLegend,'EdgeColor',[1 1 1]); %set the edge colour of the legend to white
set(gca, 'FontName', 'Helvetica')
set(gca, 'Box', 'on', 'TickDir', 'out', 'TickLength', [.02 .02], ...
    'XMinorTick', 'on', 'YMinorTick', 'on', 'YGrid', 'on', 'XGrid', 'on',...
    'XColor', [.3 .3 .3], 'YColor', [.3 .3 .3], ...
    'LineWidth', 1)
set(gcf,'color','w'); 
set(gca,'color','w'); 
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
x=1e4:1000:1e7;
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
hLegend=legend([err_ben_m],...
    'Bending test on SM45(\sigma_m=196 MPa)','location','northwest');
set(hLegend, 'FontSize', 28);
set(hLegend,'Box','on');
set(hLegend,'EdgeColor',[1 1 1]); %set the edge colour of the legend to white
set(gca, 'FontName', 'Helvetica')
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


figure(1);
saveas(gcf,'F:\Git\Anew\figures\SM45C_b1D_m_sn.png');
figure(2);
saveas(gcf,'F:\Git\Anew\figures\SM45C_b1D_m_err.png');

sp=actxserver('SAPI.SpVoice');
sp.Speak('done');
