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

NFben=[1.00E+05 2.00E+05 3.00E+05 4.00E+05 5.00E+05 6.00E+05  7.00E+05 8.00E+05 9.00E+05 1.00E+06] ;
stressben=1e6.*[311.30 289.36 276.53 267.43 260.37 254.60 249.72 245.49 241.77 238.43;...
    281.77 267.18 258.64 252.58 247.89 244.05 240.80 237.99 235.51 233.29; ...
    257.98 242.47 233.40 226.96 221.97 217.89 214.44 211.46 208.82 206.47;...
    251.82 224.42 208.40 197.03 188.21 181.00 174.91 169.63 164.97 160.81; ];%to get Smaxben

for  j=1:4 %-----------4 mean stress-----------------
    sigm=m(j);
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
            dev1=tensor-hydro(n)*eye(3)-scentre;
            dev11=dev1(1,1); dev12=dev1(1,2); dev13=dev1(1,3);
            dev21=dev1(2,1); dev22=dev1(2,2); dev23=dev1(2,3);
            dev31=dev1(3,1); dev32=dev1(3,2); dev33=dev1(3,3);
            tensor = [stressben(j,i)*sind((n+1)*360/stepnumber)+m(j,i), 0, 0 ;...
                0, 0,  0 ;...
                0, 0,  0 ; ];
            run('Damiter2.m')
            n=n+1;
        end
        Smax_ben(j,i)=max(Smax); %max sqrt of J2,a
        alp_ben(j,i)=mean(alp_ref);
        hydroplus(j,i)=1/sqrt(2)*(max(hydro)-mean(hydro))+mean(hydro);
        hydrominus(j,i)=1/sqrt(2)*(min(hydro)-mean(hydro))+mean(hydro);
    end

end

% % -------------fit directly on 1d with mean stress(1 lam)
NFben_tensor=repmat(NFben,4,1);
fun_analytical=@(parameters)parameters(1).*E.*(E+k.*nu).*parameters(2).*(parameters(2)+1).*...
    (2.*(1-alp_ben).*(E-k).*(1+nu).*(parameters(2)-1)).^-1.*...
    (Smax_ben.^(parameters(2)+1).*(y-parameters(3).*hydroplus).^(1-parameters(2))...
    +Smax_ben.^(parameters(2)+1).*(y-parameters(3).*hydrominus).^(1-parameters(2))).^-1-NFben_tensor;

% % %%-------------lsqnonlin locally stuck, issue new parameters of the non-continuous function using random
% % %%-------------generated new values-------------
% parameters=[1e5,3,0.5];
% lb  =  [1e3,   1,0];
% ub =  [1e10,   10,1];
% [parameters_fit,resnorm,exitflag]=lsqnonlin(fun_analytical,parameters,lb,ub);
% [W0,b,lamplus,lamminus]=deal(parameters_fit(1),parameters_fit(2),...
%     parameters_fit(3),parameters_fit(3));
% least_square=resnorm;
% W0_test(2)=parameters_fit(1); %original fit
% b_test(2)=parameters_fit(2); %original fit
% lamplus_test(2)=parameters_fit(3);
% lamminus_test(2)=parameters_fit(3);
% %------iterate to b converge-----------------
% p=2;
% while   resnorm>4e12
%     b=b_test(p)
%     W0=W0_test(p)
%     lamplus=lamplus_test(p)
%     lamminus=lamminus_test(p)
%     for  i=1:length(NFben)
%         n=1;
%         tensor = [stressben(i)*sind(n*360/stepnumber)+m(i) 0 0 ;...
%             0 0 0 ;...
%             0 0 0 ];
%         sigm=m(i);
%         scentre=[2*sigm/3            0                0 ;...
%             0             -sigm/3                0 ;...
%             0                       0      -sigm/3];
%         run('Damiter1.m')
%         while n<cycles*stepnumber
%             tensor = [stressben(i)*sind(n*360/stepnumber)+m(i), 0, 0 ;...
%                 0, 0, 0 ;...
%                 0, 0, 0 ; ];
%             hydro=1/3*trace(tensor);
%             dev1=tensor-hydro*eye(3)-scentre;
%             dev11=dev1(1,1); dev12=dev1(1,2); dev13=dev1(1,3);%to give \dot{dev\Sigma}dt
%             dev21=dev1(2,1); dev22=dev1(2,2); dev23=dev1(2,3);
%             dev31=dev1(3,1); dev32=dev1(3,2); dev33=dev1(3,3);
%             tensor = [stressben(i)*sind((n+1)*360/stepnumber)+m(i), 0, 0 ;...
%                 0, 0,  0 ;...
%                 0, 0,  0 ; ];
%             run('Damiter2.m')
%             n=n+1;
%         end
%         Smax_ben(i)=max(Smax);
%         alp_ben(i)=mean(alp_ref);
%         hydroplus(i)=1/sqrt(2)*(max(hydro)-mean(hydro))+mean(hydro);
%         hydrominus(i)=1/sqrt(2)*(min(hydro)-mean(hydro))+mean(hydro);
%     end
%     parameters0=[W0,b,lamplus,lamminus];  %from pth fit
%     [parameters_fit,resnorm,exitflag]=lsqnonlin(fun_analytical,parameters0,lb,ub);
%     least_square(p)=resnorm;
%     resnorm(end)
%     p=p+1;
%     W0_test(p)=parameters_fit(1);
%     b_test(p)=parameters_fit(2);
%     lamplus_test(p)=parameters_fit(3);
%     lamminus_test(p)=parameters_fit(3);
%     if abs(b_test(p)-b)/b<1e-5 %if new fit stucks
%         b_test(p)=b_test(p)*randi([110,200],1,1)/100; %--multiplied by 1.1~2
%             if b_test(p)>ub(2)
%             b_test(p)=1.01;
%         end    
%     end
%     if abs(W0_test(p)-W0)/W0<1e-5 %if new fit stucks
%         W0_test(p)=W0_test(p)*randi([110,900],1,1)/100; %--multiplied by 1.1~9
%         if W0_test(p)>ub(1)
%             W0_test(p)=lb(1);
%         end    
%     end
% end
% W0=W0_test(p);
% b=b_test(p);
% lamplus=lamplus_test(p);
% lamminus=lamminus_test(p);
% least_square=resnorm(end)
% W0
% b
% lamplus
% lamminus
% % lamminus=lamratio.*lamplus;
% save('HNAP.mat','W0','b','lamplus','lamminus','-append');
% %%%%--------------------------------------------------------------------------

parameters=[W0,b,lamplus,lamminus];
NFben_num=fun_analytical(parameters)+NFben_tensor;

for  j=1:4
    figure(1);
    hold on;
    ben_exp(j)=semilogx(NFben,Smax_ben(j,:),marker(j),'MarkerSize',12,'LineWidth', 3,'MarkerEdgeColor',color(j,:), 'MarkerFaceColor','none');
    ben_num(j)=semilogx(NFben_num(j,:),Smax_ben(j,:),linestyle(j),'color',color(j,:),'LineWidth', 3);
    figure(2);
    err_ben_m(j) = loglog(NFben,NFben_num(j,:),marker(j),'MarkerSize',12,'LineWidth', 3,'MarkerEdgeColor',color(j,:), 'MarkerFaceColor','none');
    hold on;
end

%%
%
% figure(3);
% grid on;
% set(gca ,'FontSize',30);
% hXLabel = xlabel('n','Fontsize',30, 'FontWeight' , 'bold');
% hYLabel = ylabel('S_{max}','Fontsize',30, 'FontWeight' , 'bold');
% set(gcf,'color','w'); %set figure background transparent
% set(gca,'color','w'); %set axis transparent
% set(gcf,'outerposition',get(0,'screensize'));
% set(gcf, 'PaperPositionMode', 'manual');
% set(gcf, 'PaperUnits', 'points'); %[ {inches} | centimeters | normalized | points ]
% set(gcf, 'PaperPosition', [0 0 800 800]); %set(gcf,'PaperPosition',[left,bottom,width,height])
%
% figure(4);
% grid on;
% set(gca ,'FontSize',30);
% hXLabel = xlabel('n','Fontsize',30, 'FontWeight' , 'bold');
% hYLabel = ylabel('Hydro_{+}(solid) and hydro_{-}(void)','Fontsize',30, 'FontWeight' , 'bold');
% set(gcf,'color','w'); %set figure background transparent
% set(gca,'color','w'); %set axis transparent
% set(gcf,'outerposition',get(0,'screensize'));
% set(gcf, 'PaperPositionMode', 'manual');
% set(gcf, 'PaperUnits', 'points'); %[ {inches} | centimeters | normalized | points ]
% set(gcf, 'PaperPosition', [0 0 800 800]); %set(gcf,'PaperPosition',[left,bottom,width,height])


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
% figure(3);
% saveas(gcf,'F:\Git\Anew\figures\10HNAP_b1D_m_Smax.png');
% figure(4);
% saveas(gcf,'F:\Git\Anew\figures\10HNAP_b1D_m_hydro.png');

sp=actxserver('SAPI.SpVoice');
sp.Speak('done');


