clc;
clear;
close all;
format long;
% steel SM45C data from Jabbado thesis
%-------------------------------------------------bending--------------------------------
load('SM45C.mat');
lamplus=1
lamminus=lamratio.*lamplus
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
    alp_ben_m(i)=mean(alp);
    hydroplus(i)=1/sqrt(2)*(max(hydro)-mean(hydro))+mean(hydro);
    hydrominus(i)=1/sqrt(2)*(min(hydro)-mean(hydro))+mean(hydro);
end

fun_analytical=@(parameters)parameters(1).*E.*(E+k.*nu).*parameters(2).*(parameters(2)+1).*...
    (2.*(1-alp_ben_m).*(E-k).*(1+nu).*(parameters(2)-1)).^-1.*...
    (Smax_ben_m.^(parameters(2)+1).*(y-parameters(3).*hydroplus).^(1-parameters(2))...
    +Smax_ben_m.^(parameters(2)+1).*(y-parameters(3).*hydrominus).^(1-parameters(2))).^-1-NFben';
%%------identify b with 1D_m------------
% parameters=[1,1.1,0.1];
% lb  =  [1e6,     1,0];
% ub =  [1e8,   10,1];
% [parameters_fit,resnorm,exitflag]=lsqnonlin(fun_analytical,parameters,lb,ub);
% [W0,b,lamplus,lamminus]=deal(parameters_fit(1),parameters_fit(2),parameters_fit(3),parameters_fit(3));
% least_square=resnorm;
% W0_test(2)=parameters_fit(1); %original fit
% b_test(2)=parameters_fit(2); %original fit
% lamplus_test(2)=parameters_fit(3);
% %------iterate to b converge-----------------
% p=2;
% % while abs(b_test(p)-b_test(p-1))/b_test(p-1)>1e-3
% while   resnorm>6.36e10
%     b=b_test(p)
%     W0=W0_test(p)
%     lamplus=lamplus_test(p)
%     % save('SM45C.mat','W0','b','-append');
%     for  i=1:length(NFben)
%         n=1;
%         tensor = [stressben(i)*sind(n*360/stepnumber)+m 0 0 ;...
%             0 0 0 ;...
%             0 0 0 ];
%         sigm=m;
%         scentre=[2*sigm/3            0                0 ;...
%             0             -sigm/3                0 ;...
%             0                       0      -sigm/3];
%         run('Damiter1.m')
%         while n<cycles*stepnumber
%             tensor = [stressben(i)*sind(n*360/stepnumber)+m, 0, 0 ;...
%                 0, 0, 0 ;...
%                 0, 0, 0 ; ];
%             hydro=1/3*trace(tensor);
%             dev1=tensor-hydro*eye(3)-scentre;
%             dev11=dev1(1,1); dev12=dev1(1,2); dev13=dev1(1,3);%to give \dot{dev\Sigma}dt
%             dev21=dev1(2,1); dev22=dev1(2,2); dev23=dev1(2,3);
%             dev31=dev1(3,1); dev32=dev1(3,2); dev33=dev1(3,3);
%             tensor = [stressben(i)*sind((n+1)*360/stepnumber)+m, 0, 0 ;...
%                 0, 0,  0 ;...
%                 0, 0,  0 ; ];
%             run('Damiter2.m')
%             n=n+1;
%         end
%         Smax_ben_m(i)=max(Smax);
%         alp_ben_m(i)=mean(alp);
%         hydroplus(i)=1/sqrt(2)*(max(hydro)-mean(hydro))+mean(hydro);
%         hydrominus(i)=1/sqrt(2)*(min(hydro)-mean(hydro))+mean(hydro);
%     end
%     parameters0=[W0,b,lamplus];  %from pth fit
%     [parameters_fit,resnorm,exitflag]=lsqnonlin(fun_analytical,parameters0,lb,ub);
%     least_square(p)=resnorm;
%     resnorm(end)
%     p=p+1;
%     W0_test(p)=parameters_fit(1);
%     b_test(p)=parameters_fit(2);
%     lamplus_test(p)=parameters_fit(3);
%     if abs(b_test(p)-b)/b<1e-4 %if new fit stucks
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
% 
% W0=W0_test(p);
% b=b_test(p);
% lamplus=lamplus_test(p);
% W0
% b
% lamplus
% lamminus=lamratio.*lamplus;
% 
% save('SM45C.mat','lamplus','lamminus','W0','b','-append');
parameters=[W0,b,lamplus];
NF_num=fun_analytical(parameters)+NFben';

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

save('SM45C.mat','lamplus','lamminus','-append');

figure(1);
saveas(gcf,'F:\Git\Anew\figures\SM45C_b1D_m_sn.png');
figure(2);
saveas(gcf,'F:\Git\Anew\figures\SM45C_b1D_m_err.png');

sp=actxserver('SAPI.SpVoice');
sp.Speak('done');
