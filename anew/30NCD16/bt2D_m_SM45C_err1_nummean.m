clc;
clear;
close all;
% steel SM45C data from Jabbado thesis
%-------------------------------------------------bending--------------------------------
y=638e6;
k=1e9;
E=213e9;
nu=0.3;                     %poisson's ratio
a=0.01;
W0=13e6;
lam=0.1;
b=1.2;
fb=b;
m=196e6;%mean tension
stepnumber=100;        %devide one cycle in 100 parts
x= [0.999305042	0.996340117	0.991013371	0.983336254	0.973326828	0.9610088	0.946411375	0.929569172	0.910522137...
    0.889315446	0.865999398	0.840629296	0.813265315	0.783972359	0.752819907	0.71988185	0.685236313	0.648965471...
    0.611155355	0.571895646	0.531279464	0.489403146	0.446366017	0.402270158	0.357220158	0.311322872	0.264687162...
    0.217423644	0.16964442	0.121462819	0.072993122	0.024350293	-0.024350293	-0.072993122	-0.121462819	-0.16964442...
    -0.217423644	-0.264687162	-0.311322872	-0.357220158	-0.402270158	-0.446366017	-0.489403146	-0.531279464...
    -0.571895646	-0.611155355	-0.648965471	-0.685236313	-0.71988185	-0.752819907	-0.783972359	-0.813265315...
    -0.840629296	-0.865999398	-0.889315446	-0.910522137	-0.929569172	-0.946411375	-0.9610088	-0.973326828...
    -0.983336254	-0.991013371	-0.996340117	-0.999305042];
weight=[0.001783281	0.004147033	0.006504458	0.00884676	0.011168139	0.013463048	0.01572603	0.017951716	0.020134823...
    0.022270174	0.024352703	0.02637747	0.028339673	0.030234657	0.032057928	0.033805162	0.035472213	0.037055129	0.038550153...
    0.039953741	0.041262563	0.042473515	0.043583725	0.044590558	0.045491628	0.046284797	0.046968183	0.047540166	0.047999389...
    0.048344762	0.048575467	0.048690957	0.048690957	0.048575467	0.048344762	0.047999389	0.047540166	0.046968183	0.046284797...
    0.045491628	0.044590558	0.043583725	0.042473515	0.041262563	0.039953741	0.038550153	0.037055129	0.035472213	0.033805162...
    0.032057928	0.030234657	0.028339673	0.02637747	0.024352703	0.022270174	0.020134823	0.017951716	0.01572603	0.013463048...
    0.011168139	0.00884676	0.006504458	0.004147033	0.001783281];
%---------------a b correspond to high and low bounds of cyclic load, 1 2
%correspond to  different parameters comparison------------------
NF=[53E3 59.2E3 70.1E3 86.3E3 89.9E3 92.1E3 102E3 135E3 351E3 394E3] ;
stressben=[441E6 286E6 464E6 473E6 173E6 403E6 437E6 167E6 357E6 182E6];%to get Smaxben
stresstor=[215E6 309E6 155E6 136E6 334E6 209E6 177E6 321E6 179E6 274E6];%to get Smaxtor
%---------------------Numerical to get the mean value via several cycles-----------------------------
tic;
% for  i=1:length(NF)
for  i=1:length(NF)
    D=1e-16;
    n=1;       %initial recording point
    tensor = [stressben(i)*sind(n*360/stepnumber)+m -stresstor(i)*sind(n*360/stepnumber) 0 ;...
        -stresstor(i)*sind(n*360/stepnumber) 0  0 ;...
        0 0 0 ];
    %---------------------to get the the first Sb-----------------------------
    hydro=1/3*trace(tensor);
    yield=y-lam*hydro; %micro yield strength at n=1
    dev1=tensor-hydro*eye(3);
    dev11=dev1(1,1); dev12=dev1(1,2); dev13=dev1(1,3);
    dev21=dev1(2,1); dev22=dev1(2,2); dev23=dev1(2,3);
    dev31=dev1(3,1); dev32=dev1(3,2); dev33=dev1(3,3);
    Smax=norm(dev1,'fro');
    s= (x/2+1/2).^(1/(1-b)); %1*64
    
    trial11=dev11; trial12=dev12; trial13=dev13;
    trial21=dev21; trial22=dev22; trial23=dev23;
    trial31=dev31; trial32=dev32; trial33=dev33;
    
    normtrial(1)=norm([trial11, trial12, trial13; trial21, trial22, trial23;trial31, trial32, trial33],'fro');
    
    eta=bsxfun(@minus,bsxfun(@times,normtrial(1)/yield,s),1); %compare normtrial with yield/s
    eta(eta<0)=0; %only keep normtrials which are larger than yield/s
    
    Sb11=bsxfun(@rdivide,trial11,bsxfun(@plus,eta,1));Sb12=bsxfun(@rdivide,trial12,bsxfun(@plus,eta,1));Sb13=bsxfun(@rdivide,trial13,bsxfun(@plus,eta,1));
    Sb21=bsxfun(@rdivide,trial21,bsxfun(@plus,eta,1));Sb22=bsxfun(@rdivide,trial22,bsxfun(@plus,eta,1));Sb23=bsxfun(@rdivide,trial23,bsxfun(@plus,eta,1));
    Sb31=bsxfun(@rdivide,trial31,bsxfun(@plus,eta,1));Sb32=bsxfun(@rdivide,trial32,bsxfun(@plus,eta,1));Sb33=bsxfun(@rdivide,trial33,bsxfun(@plus,eta,1));
    %1*64 for each Sb element
    Sbtensor=[Sb11; Sb12; Sb13; Sb21; Sb22; Sb23;Sb31; Sb32; Sb33];
    normSb=sqrt(sum(Sbtensor.^2)); %sum(a) sums all the colume
    
    % existsOnGPU(normSb)
    Ws=(bsxfun(@minus,normtrial,bsxfun(@rdivide, yield,s))<=0).*...
        (0)+...
        (bsxfun(@minus,normtrial,bsxfun(@rdivide, yield,s))>0).*...
        ((E-k)*(1+nu)*(2*E*(E+k*nu))^-1*bsxfun(@times,weight,bsxfun(@rdivide,bsxfun(@times,bsxfun(@minus,normtrial,bsxfun(@rdivide, yield,s)),yield),s)));
    
    W(1)= sum(Ws);
    sequence=((Smax*yield^-1)*(1-Smax*yield^-1)^-1)^b;
    sequence(sequence<0)=0;
    alp(1)=1-a*sequence;
    D=D+D^alp(1)*W(1)/W0;
    
    while n<1*stepnumber
        tensor = [stressben(i)*sind(n*360/stepnumber)+m, -stresstor(i)*sind(n*360/stepnumber), 0 ;...
            -stresstor(i)*sind(n*360/stepnumber), 0,  0 ;...
            0,0, 0; ];
        hydro=1/3*trace(tensor);
        dev1=tensor-hydro*eye(3);
        dev11=dev1(1,1); dev12=dev1(1,2); dev13=dev1(1,3);%to give \dot{dev\Sigma}dt
        dev21=dev1(2,1); dev22=dev1(2,2); dev23=dev1(2,3);
        dev31=dev1(3,1); dev32=dev1(3,2); dev33=dev1(3,3);
        
        tensor = [stressben(i)*sind((n+1)*360/stepnumber)+m, -stresstor(i)*sind((n+1)*360/stepnumber), 0 ;...
            -stresstor(i)*sind((n+1)*360/stepnumber), 0,  0 ;...
            0, 0, 0; ];
        hydro=1/3*trace(tensor);
        yield(n+1)=y-lam*hydro; %yield stress at time step n+1
        devn=tensor-hydro*eye(3);
        dev11g=devn(1,1); dev12g=devn(1,2); dev13g=devn(1,3);
        dev21g=devn(2,1); dev22g=devn(2,2); dev23g=devn(2,3);
        dev31g=devn(3,1); dev32g=devn(3,2); dev33g=devn(3,3);
        Smax(n+1)=norm(devn,'fro');
        
        trial11=bsxfun(@plus,Sb11,(dev11g-dev11)); trial12=bsxfun(@plus,Sb12,(dev12g-dev12));trial13=bsxfun(@plus,Sb13,(dev13g-dev13));
        trial21=bsxfun(@plus,Sb21,(dev21g-dev21)); trial22=bsxfun(@plus,Sb22,(dev22g-dev22));trial23=bsxfun(@plus,Sb23,(dev23g-dev23));
        trial31=bsxfun(@plus,Sb31,(dev31g-dev31)); trial32=bsxfun(@plus,Sb32,(dev32g-dev32));trial33=bsxfun(@plus,Sb33,(dev33g-dev33));
        trialtensor=[trial11; trial12; trial13; trial21; trial22; trial23;trial31; trial32; trial33];
        normtrial=sqrt(sum(trialtensor.^2));
        
        eta=bsxfun(@minus,bsxfun(@times,normtrial/yield(n+1),s),1); %1*64
        eta(eta<0)=0;
        
        Sb11=bsxfun(@rdivide,trial11,bsxfun(@plus,eta,1));Sb12=bsxfun(@rdivide,trial12,bsxfun(@plus,eta,1));Sb13=bsxfun(@rdivide,trial13,bsxfun(@plus,eta,1));
        Sb21=bsxfun(@rdivide,trial21,bsxfun(@plus,eta,1));Sb22=bsxfun(@rdivide,trial22,bsxfun(@plus,eta,1));Sb23=bsxfun(@rdivide,trial23,bsxfun(@plus,eta,1));
        Sb31=bsxfun(@rdivide,trial31,bsxfun(@plus,eta,1));Sb32=bsxfun(@rdivide,trial32,bsxfun(@plus,eta,1));Sb33=bsxfun(@rdivide,trial33,bsxfun(@plus,eta,1));
        %1*64 for each Sb element
        Sbtensor=[Sb11; Sb12; Sb13; Sb21; Sb22; Sb23;Sb31; Sb32; Sb33];
        normSb=sqrt(sum(Sbtensor.^2)); %sum(a) sums all the colume
        
        Ws=(bsxfun(@minus,normtrial,bsxfun(@rdivide, yield(n+1),s))<=0).*...
            (0)+...
            (bsxfun(@minus,normtrial,bsxfun(@rdivide, yield(n+1),s))>0).*...
            ((E-k)*(1+nu)*(2*E*(E+k*nu))^-1*bsxfun(@times,weight,bsxfun(@rdivide,bsxfun(@times,bsxfun(@minus,normtrial,bsxfun(@rdivide, yield(n+1),s)),yield(n+1)),s)));
        W(n+1)= sum(Ws);
        
        sequence=((Smax(n+1)*yield(n+1)^-1)*(1-Smax(n+1)*yield(n+1)^-1)^-1)^b;
        sequence(sequence<0)=0;
        alp(n+1)=1-a*sequence;
        D=D+D^alp(n+1)*W(n+1)/W0;
        
%         figure(3);
%         hold on;
%         yield1=plot (n+1,yield(n+1)*s(20).^-1, 'LineStyle', 'none','LineWidth', 1, 'Marker', 'o', 'MarkerSize', 8, ...
%             'MarkerEdgeColor',  'none', 'MarkerFaceColor' , 'c');
%         Trial1=plot (n+1,normtrial(20),'LineStyle', 'none','LineWidth', 1,'Marker', '^', 'MarkerSize', 10, ...
%             'MarkerEdgeColor','r', 'MarkerFaceColor','r');
%         Sb1=plot (n+1,normSb(20),'LineStyle', 'none','LineWidth', 1,'Marker', 'v', 'MarkerSize', 10, ...
%             'MarkerEdgeColor','g', 'MarkerFaceColor','m');
%         dev=plot (n+1,Smax(n+1),'LineStyle', 'none','LineWidth', 1,'Marker', 's', 'MarkerSize', 11, ...
%             'MarkerEdgeColor','none', 'MarkerFaceColor',[238 18 137]/255);
%         yield8=plot (n+1,yield(n+1)*s(12).^-1,'LineStyle', 'none','LineWidth', 1,'Marker', 'o', 'MarkerSize', 8, ...
%             'MarkerEdgeColor', 'none', 'MarkerFaceColor', 'b');
%         Trial8=plot (n+1,normtrial(12),'LineStyle', 'none','LineWidth', 1,'Marker', '^', 'MarkerSize', 10, ...
%             'MarkerEdgeColor', [1 0.5 0], 'MarkerFaceColor',[1 0.5 0]);
%         Sb8=plot (n+1,normSb(12),'LineStyle', 'none','LineWidth', 1,'Marker', 'v', 'MarkerSize', 10, ...
%             'MarkerEdgeColor','k', 'MarkerFaceColor','k');
%         figure(4);
%         hold on;
%         plot_alpha=plot (n+1,alp(n+1),'LineStyle', 'none','LineWidth', 1,'Marker', 's', 'MarkerSize', 10, ...
%             'MarkerEdgeColor','k', 'MarkerFaceColor','k');
%         figure(5);
%         hold on;
%         plot_W=plot (n+1,W(n+1),'LineStyle', 'none','LineWidth', 1,'Marker', 'o', 'MarkerSize', 10, ...
%             'MarkerEdgeColor','m', 'MarkerFaceColor','m');
        n=n+1;
    end
    %control lam to make min(sminmin)>1
    sminmin(i)=min(yield(1:n)/Smax(1:n));
    Smax_m(i)=mean(Smax);
    Wcyc_m(i)=mean(W);
    alpbt_m(i)=mean(alp);
    NF_num(i)=W0*((1-alpbt_m(i))*Wcyc_m(i)).^-1;
end
toc;
% figure(3);
% set(gca ,'FontSize',30);
% hXLabel = xlabel('NF_{num}','Fontsize',30, 'FontWeight' , 'bold');
% hYLabel = ylabel('Stress','Fontsize',30, 'FontWeight' , 'bold');
% set(gcf,'outerposition',get(0,'screensize'));
% figure(4);
% set(gca ,'FontSize',30);
% hXLabel = xlabel('NF_{num}','Fontsize',30, 'FontWeight' , 'bold');
% hYLabel = ylabel('\alpha','Fontsize',30, 'FontWeight' , 'bold');
% set(gcf,'outerposition',get(0,'screensize'));
% figure(5);
% set(gca ,'FontSize',30);
% hXLabel = xlabel('NF_{num}','Fontsize',30, 'FontWeight' , 'bold');
% hYLabel = ylabel('W','Fontsize',30, 'FontWeight' , 'bold');
% set(gcf,'outerposition',get(0,'screensize'));
sminmin
alpbt_m
%% 
%--------use average smin to get NF--------
figure(1);
smax_exp=semilogx(NF(1:i),Smax_m(1:i),'o','MarkerSize',12,'LineWidth', 3,'MarkerEdgeColor',[208 32 144]/255, 'MarkerFaceColor','none');
hold on;
smax_num=semilogx(NF_num(1:i),Smax_m(1:i),'-r','LineWidth', 3);

set(gca ,'FontSize',30);
hXLabel = xlabel('NF_{num}','Fontsize',30, 'FontWeight' , 'bold');
hYLabel = ylabel('S_{max}','Fontsize',30, 'FontWeight' , 'bold');
hLegend=legend([smax_exp,smax_num],...
    'Bending-torsion experimental result','Bending-torsion numerical result','location','bestoutside');
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
set(gcf, 'PaperPosition', [0 0 1800 1000]); %set(gcf,'PaperPosition',[left,bottom,width,height])

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
    'Bending-torsion 90 degree out-of-phase test with mean stress on SM45C','location','bestoutside');
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
set(gcf, 'PaperPosition', [0 0 1800 1000]); %set(gcf,'PaperPosition',[left,bottom,width,height])
% figure(1);
% saveas(gcf,'F:\Git\Anew\figures\bt2D_m_SM45C_sn.png');
% figure(2);
% saveas(gcf,'F:\Git\Anew\figures\bt2D_m_SM45C_err1.png');

sp=actxserver('SAPI.SpVoice');
sp.Speak('done done done');

