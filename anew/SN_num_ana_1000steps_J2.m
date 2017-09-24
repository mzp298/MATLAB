clear;clc;close all;format long;
load('gaussian.mat');
%---------------------Verified parameters in random loading case-----------------------------
b=1.1;                    %weakening scales distribution exponent
W0=3.27e5;            %dissipated energy to failure per unit volume
%---------------------Verified parameters in constant loading case-----------------------------
y=1080e6;           %macroscopic yield stress
E=191e9;              %Young's modulus
k=1e9;                 %hardening parameter
nu=0.38;                     %poisson's ratio
lam=0.3;               %hydrostatic pressure sensitivity
fb=b;                %Major damage power
%---------------------Verified parameters in constant loading case-----------------------------
stepnumber=1000;        %devide one cycle in 200 parts
delta_alp=0.01;       %optimal alpha filter threshold
cycles=2;          %numerical cycles to adaptation
m=0;                   % mean stress
load = [y/1.2 y/1.3 y/1.5 y/2 y/2.5 y/3 y/3.5 y/4 y/4.5 y/5 y/5.5 y/6];
% hydro=1/3*sum(max(load)+0+0);
% yield=y-lam*hydro;
% a=2.5*(yield/(sqrt(1/3)*max(load))-1)^b
a=1e-3

% load=1000:1000:y; %---------to plot the SN curve line------------

%---------------------1 Numerical method with optimal time steps-----------------------------
for  i=1:length(load) %experimental points index
    n=1;
    tensor = [load(i)*sind(n*360/stepnumber)+m 0 0 ;...
        0 0 0 ;...
        0 0 0 ];
    run('F:\Git\MATLAB\anew\Damiter1.m')
    g=1; %first reference alp, W, n index when generate
    while n<cycles*stepnumber
        tensor = [load(i)*sind(n*360/stepnumber)+m, 0, 0 ;...
            0, 0, 0 ;...
            0, 0, 0 ; ];
        hydro=1/3*trace(tensor);
        dev1=tensor-hydro*eye(3);
        dev11=dev1(1,1); dev12=dev1(1,2); dev13=dev1(1,3);
        dev21=dev1(2,1); dev22=dev1(2,2); dev23=dev1(2,3);
        dev31=dev1(3,1); dev32=dev1(3,2); dev33=dev1(3,3);
        tensor = [load(i)*sind((n+1)*360/stepnumber)+m, 0, 0 ;...
            0, 0,  0 ;...
            0, 0,  0 ; ];
        run('F:\Git\MATLAB\anew\Damiter2.m')
        n=n+1;
    end
    alp_m(i)=mean(alp)
%     figure(3);
%     plot(alp,'b','Marker','s','MarkerSize',20,'MarkerFaceColor','none','LineStyle','none', 'LineWidth', 2);
    %----------scalar iteration
    j=ceil(length(alp_ref)/cycles); %first reference alp, W, n index when iterate(after adaptation)
    D=1e-16;
    n=1; %first D index
    while D<1 %-----------the optimal time steps can be iterated with scalar
        D=D+D^alp_ref(j)*W_ref(j)/W0; % explicit(require D \neq 0, normally D=1e-16)
        j=j+1;
        if j+1>=length(alp_ref)
            j=ceil(length(alp_ref)/cycles);
        end
        n=n+1;
    end
    NF_num_opt(i)=(n-1)/(length(n_ref)/cycles);
end
reduction_factor=length(alp_ref)/cycles/stepnumber % to see how much optimal timesteps comparing to the original
%
%% 
% 
%---------------------2 Numerical method-----------------------------
for i=1:length(load)
    D=1e-16;
% D=0;
    n=1;       %initial recording point
    %---------------------to get the the first Sb-----------------------------
    loadtensor= [load(i) 0 0;0 0 0;0 0 0];
    stress11=m+load(i)*sind(n*360/stepnumber);
    hydro=1/3*sum(stress11(1)+0+0);
    yield=y-lam*hydro; %macro yield strength considering mean stress effect
    dev1=[stress11 0 0 ;0 0 0 ;0 0 0 ]-hydro*eye(3);
    dev11=dev1(1,1); dev12=dev1(1,3); dev13=dev1(1,3);
    dev21=dev1(2,1); dev22=dev1(2,2); dev23=dev1(2,3);
    dev31=dev1(3,1); dev32=dev1(3,2); dev33=dev1(3,3);
    Smax=1/sqrt(2).*norm(dev1,'fro');
    s= (x/2+1/2).^(1/(1-b)); %1*64
    
    trial11=dev11; trial12=dev12; trial13=dev13;
    trial21=dev21; trial22=dev22; trial23=dev23;
    trial31=dev31; trial32=dev32; trial33=dev33;
    
    Smaxtrial(1)=1/sqrt(2).*norm([trial11, trial12, trial13; trial21, trial22, trial23;trial31, trial32, trial33],'fro');
    
    eta=bsxfun(@minus,bsxfun(@times,Smaxtrial(1)/yield,s),1); %compare Smaxtrial with yield/s
    eta(eta<0)=0; %only keep Smaxtrials which are larger than yield/s
    
    Sb11=bsxfun(@rdivide,trial11,bsxfun(@plus,eta,1));Sb12=bsxfun(@rdivide,trial12,bsxfun(@plus,eta,1));Sb13=bsxfun(@rdivide,trial13,bsxfun(@plus,eta,1));
    Sb21=bsxfun(@rdivide,trial21,bsxfun(@plus,eta,1));Sb22=bsxfun(@rdivide,trial22,bsxfun(@plus,eta,1));Sb23=bsxfun(@rdivide,trial23,bsxfun(@plus,eta,1));
    Sb31=bsxfun(@rdivide,trial31,bsxfun(@plus,eta,1));Sb32=bsxfun(@rdivide,trial32,bsxfun(@plus,eta,1));Sb33=bsxfun(@rdivide,trial33,bsxfun(@plus,eta,1));

    Ws=(bsxfun(@minus,Smaxtrial,bsxfun(@rdivide, yield,s))<=0).*...
        (0)+...
        (bsxfun(@minus,Smaxtrial,bsxfun(@rdivide, yield,s))>0).*...
        ((E-k)*(1+nu)*(2*E*(E+k*nu))^-1*bsxfun(@times,weight,bsxfun(@rdivide,bsxfun(@times,bsxfun(@minus,Smaxtrial,bsxfun(@rdivide, yield,s)),yield),s)));
    
    W= sum(Ws);
    sequence=((Smax*yield^-1)*(1-Smax*yield^-1)^-1)^b;
    
    sequence(sequence<0)=0;
    alp=1-a*sequence;
    D=D+D^alp*W/W0;
    

    
% D= (D^(1-alp)+(1-alp).*W/W0).^(1/(1-alp)); % implicit
    while  D<1 %-----------the optimal time steps can be iterated with scalar
        stress11=m+load(i)*sind(n*360/stepnumber);
        hydro=1/3*sum(stress11+0+0);
        yield=y-lam*hydro;
        dev1=[stress11 0 0 ;0 0 0 ;0 0 0 ]-hydro*eye(3);
        dev11=dev1(1,1); dev12=dev1(1,2); dev13=dev1(1,3);
        dev21=dev1(2,1); dev22=dev1(2,2); dev23=dev1(2,3);
        dev31=dev1(3,1); dev32=dev1(3,2); dev33=dev1(3,3);
        
        stress11=m+load(i)*sind((n+1)*360/stepnumber);
        hydro=1/3*sum(stress11+0+0);
        devn=[stress11 0 0;0 0 0;0 0 0]-hydro*eye(3);
        dev11g=devn(1,1); dev12g=devn(1,2); dev13g=devn(1,3);
        dev21g=devn(2,1); dev22g=devn(2,2); dev23g=devn(2,3);
        dev31g=devn(3,1); dev32g=devn(3,2); dev33g=devn(3,3);
        Smax=1/sqrt(2).*norm(devn,'fro');
        
        trial11=bsxfun(@plus,Sb11,(dev11g-dev11)); trial12=bsxfun(@plus,Sb12,(dev12g-dev12));trial13=bsxfun(@plus,Sb13,(dev13g-dev13));
        trial21=bsxfun(@plus,Sb21,(dev21g-dev21)); trial22=bsxfun(@plus,Sb22,(dev22g-dev22));trial23=bsxfun(@plus,Sb23,(dev23g-dev23));
        trial31=bsxfun(@plus,Sb31,(dev31g-dev31)); trial32=bsxfun(@plus,Sb32,(dev32g-dev32));trial33=bsxfun(@plus,Sb33,(dev33g-dev33));
        trialtensor=[trial11; trial12; trial13; trial21; trial22; trial23;trial31; trial32; trial33];
        Smaxtrial=1/sqrt(2).*sqrt(sum(trialtensor.^2));
        
        eta=bsxfun(@minus,bsxfun(@times,Smaxtrial/yield,s),1); %1*64
        eta(eta<0)=0;
        
        Sb11=bsxfun(@rdivide,trial11,bsxfun(@plus,eta,1));Sb12=bsxfun(@rdivide,trial12,bsxfun(@plus,eta,1));Sb13=bsxfun(@rdivide,trial13,bsxfun(@plus,eta,1));
        Sb21=bsxfun(@rdivide,trial21,bsxfun(@plus,eta,1));Sb22=bsxfun(@rdivide,trial22,bsxfun(@plus,eta,1));Sb23=bsxfun(@rdivide,trial23,bsxfun(@plus,eta,1));
        Sb31=bsxfun(@rdivide,trial31,bsxfun(@plus,eta,1));Sb32=bsxfun(@rdivide,trial32,bsxfun(@plus,eta,1));Sb33=bsxfun(@rdivide,trial33,bsxfun(@plus,eta,1));
        
        Ws=(bsxfun(@minus,Smaxtrial,bsxfun(@rdivide, yield,s))<=0).*...
            (0)+...
            (bsxfun(@minus,Smaxtrial,bsxfun(@rdivide, yield,s))>0).*...
            ((E-k)*(1+nu)*(2*E*(E+k*nu))^-1*bsxfun(@times,weight,bsxfun(@rdivide,bsxfun(@times,bsxfun(@minus,Smaxtrial,bsxfun(@rdivide, yield,s)),yield),s)));
        
        W= sum(Ws);
        sequence=((Smax*yield^-1)*(1-Smax*yield^-1)^-1)^b;
        
        sequence(sequence<0)=0;
        alp(n+1)=1-a*sequence;
        D=D+D^alp(n+1)*W/W0;
% D= (D^(1-alp(n+1))+(1-alp(n+1)).*W/W0).^(1/(1-alp(n+1))); % implicit
        n=n+1;
    end
    NF_num(i)=n/stepnumber
end

%%
%---------------------3 Analytical calculation after integrated of D with mean alp from numerical original-----------------------------
m=0;
D0=1e-16;
for i=1:length(load)
    tensor= [load(i) 0 0;0 0 0;0 0 0];
    hydro=1/3*trace(tensor);
    yield=y-lam*hydro; %macro yield strength considering mean stress effect
    dev1=tensor-hydro*eye(3);
    Smax1(i)=1/sqrt(2).*norm(dev1,'fro');
    Wcyc1(i)=4*(E-k)*(1+nu)*(b-1)/(E*(E+k*nu)*b*(b+1)).*Smax1(i).^(b+1).*yield.^(1-b) ;
    NF_direct(i)=W0.*(Wcyc1(i).*(1-alp_m(i))).^-1*(1-D0^(1-alp_m(i)));
end
meanalp=alp_m'
optimalNF=NF_num_opt'
numericalNF=NF_num'
analyticalNF=NF_direct'
error=((NF_num(1:length(load))-NF_direct(1:length(load))).*NF_direct(1:length(load)).^-1)'

figure(1)
numerical=semilogx(NF_num,load, 's','MarkerSize',35, 'MarkerFaceColor','none', 'LineWidth', 5,'MarkerEdgeColor','b');
hold on;
%  numerical_opt=semilogx(NF_num_opt,load, 'o','MarkerSize',35, 'MarkerFaceColor','none', 'LineWidth', 5,'MarkerEdgeColor','g');
analytical_integration=semilogx(NF_direct,load,  '^','MarkerSize',35, 'MarkerFaceColor','none', 'LineWidth', 5,'MarkerEdgeColor','r');
hLegend=legend([numerical,...
%       numerical_opt,...
    analytical_integration,...
    ],...
    '1 Numerical results with original time steps ($\delta D=D^\alpha\frac{\dot{W}}{W_0}\delta t$)',...
    '2 Analytical results with $\alpha_m$ from numerical calculation ($N_F=\frac{W_0}{( 1-\alpha_m)W_{cyc}}$)',...
    'Location','best');
%     '2 Numerical results with optimal time steps ($\delta D=D^{\alpha_{ref}}\frac{W_{ref}}{W_0}\delta t_{ref}$)',...
set(hLegend,'Interpreter','latex');
set([hLegend, gca], 'FontSize', 35)
set(hLegend,'Box','on');
set(hLegend,'EdgeColor',[1 1 1]); %set the edge colour of the legend to white
grid on;
grid minor;
ylim([0 y]);
set(gca,'XTick',[1e0,1e1,1e2,1e3,1e4,1e5,1e6,1e7,1e8,1e9],'FontSize',40 );
hXLabel = xlabel('N_F');
hYLabel =ylabel('Stress');
set(gca, 'FontName', 'Helvetica')
set([hXLabel, hYLabel], 'FontName', 'AvantGarde')
set([hXLabel, hYLabel], 'FontSize',40)
set(gcf,'color','w'); %set figure background transparent
set(gca,'color','w'); %set axis transparent
set(gcf,'outerposition',get(0,'screensize'));
set(gcf, 'PaperPositionMode', 'manual');
set(gcf, 'PaperUnits', 'points'); %[ {inches} | centimeters | normalized | points ]
set(gcf, 'PaperPosition', [0 0 1600 1000]); %set(gcf,'PaperPosition',[left,bottom,width,height])

figure(2)
errorbar_plot=semilogx(NF_num,error,'k','Marker','*','MarkerSize',35,'MarkerFaceColor','none','LineStyle','none', 'LineWidth', 5);
hLegend=legend([errorbar_plot],...
    ['Relative error between analytical and numerical results',sprintf('\n'),'(a=0.001)'],...
    'Location','Southeast');
set(hLegend,'Interpreter','latex');
set([hLegend, gca], 'FontSize', 35)
set(hLegend,'Box','on');
set(hLegend,'EdgeColor',[1 1 1]); %set the edge colour of the legend to white
grid on;
grid minor;
ylim([min(error) max(error)]);
set(gca,'XTick',[1e0,1e1,1e2,1e3,1e4,1e5,1e6,1e7,1e8,1e9],'FontSize',40 );
hXLabel = xlabel('N_F');
hYLabel =ylabel('Relative error');
set(gca, 'FontName', 'Helvetica')
set([hXLabel, hYLabel], 'FontName', 'AvantGarde')
set([hXLabel, hYLabel], 'FontSize',40)
set(gcf,'color','w'); %set figure background transparent
set(gca,'color','w'); %set axis transparent
set(gcf,'outerposition',get(0,'screensize'));
set(gcf, 'PaperPositionMode', 'manual');
set(gcf, 'PaperUnits', 'points'); %[ {inches} | centimeters | normalized | points ]
set(gcf, 'PaperPosition', [0 0 1600 1000]); %set(gcf,'PaperPosition',[left,bottom,width,height])

figure(1)
saveas(gcf,'F:\Git\Anew\figures\SN_num_ana_stepnumber=1000.png');
figure(2)
saveas(gcf,'F:\Git\Anew\figures\SN_num_ana_stepnumber=1000_err.png');
