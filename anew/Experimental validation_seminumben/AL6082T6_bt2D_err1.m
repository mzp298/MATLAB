clc;
clear;
close all;
% steel AL6082T6 data from Jabbado thesis
%-------------------------------------------------bending--------------------------------
format long;
load('AL6082T6.mat'); %final fitting
load('gaussian.mat');
% vpa(x,289);
% vpa(weight,289);

NF=[266435
1119254
1416225
83000
75695
630325
157210
126470
35450
68465
252658
316149
] ;
stressben=1e6.*[57
51
51
71
70
59
132
132
144
149
122
116
];%to get Smaxben
stresstor=1e6.*[100
84
84
118
118
99
97
99
107
105
90
83
];%to get Smaxtor

m=zeros(size(NF));%mean tension
    sigm=m(1);
%---------------------Numerical to get the mean value via several cycles-----------------------------
for  i=1:length(NF)
    n=1;       %initial recording point
    tensor = [stressben(i)*sind(n*360/stepnumber)+m(i) stresstor(i)*sind(n*360/stepnumber) 0 ;...
        stresstor(i)*sind(n*360/stepnumber) 0  0 ;...
        0 0 0 ];
    %---------------------to get the the first Sb-----------------------------
    run('Damiter1.m')
    while n<cycles*stepnumber
        tensor = [stressben(i)*sind(n*360/stepnumber)+m(i), stresstor(i)*sind(n*360/stepnumber), 0 ;...
            stresstor(i)*sind(n*360/stepnumber), 0,  0 ;...
            0,0, 0; ];
        hydro(n)=1/3*trace(tensor);
        dev1=tensor-hydro(n)*eye(3)-scentre;
        dev11=dev1(1,1); dev12=dev1(1,2); dev13=dev1(1,3);%to give \dot{dev\Sigma}dt
        dev21=dev1(2,1); dev22=dev1(2,2); dev23=dev1(2,3);
        dev31=dev1(3,1); dev32=dev1(3,2); dev33=dev1(3,3);
        tensor = [stressben(i)*sind((n+1)*360/stepnumber)+m(i), stresstor(i)*sind((n+1)*360/stepnumber), 0 ;...
            stresstor(i)*sind((n+1)*360/stepnumber), 0,  0 ;...
            0, 0, 0; ];
        run('Damiter2.m')
        n=n+1;
    end
    Smax_bt2d(i)=max(Smax); %max sqrt of J2,a
    alp_bt2d(i)=mean(alp_ref);
    hydroplus(i)=1/sqrt(2)*(max(hydro)-mean(hydro))+mean(hydro);
    hydrominus(i)=1/sqrt(2)*(min(hydro)-mean(hydro))+mean(hydro);
end
fun_analytical=@(parameters)parameters(1).*E.*(E+k.*nu).*parameters(2).*(parameters(2)+1).*...
    (2.*(1-alp_bt2d).*(E-k).*(1+nu).*(parameters(2)-1)).^-1.*...
    (Smax_bt2d.^(parameters(2)+1).*(y-parameters(3).*hydroplus.*(hydroplus>0)...
    -parameters(3).*hydroplus.*(hydroplus<0)).^(1-parameters(2))...
    +Smax_bt2d.^(parameters(2)+1).*(y-parameters(3).*hydrominus.*(hydrominus<0)...
    -parameters(3).*hydrominus.*(hydrominus>0)).^(1-parameters(2))).^-1-NF';
% %%-------------lsqnonlin locally stuck, issue new parameters of the non-continuous function using random
% %%-------------generated new values-------------
parameters=[1e5,1.1,0.5];
lb  =  [1e3,   1,0];
ub =  [1e10,   10,1];
[parameters_fit,resnorm,exitflag]=lsqnonlin(fun_analytical,parameters,lb,ub);
[W0,b,lamplus,lamminus]=deal(parameters_fit(1),parameters_fit(2),parameters_fit(3),parameters_fit(3));
least_square=resnorm;
W0_test(2)=parameters_fit(1); %original fit
b_test(2)=parameters_fit(2); %original fit
lamplus_test(2)=parameters_fit(3);
lamminus_test(2)=parameters_fit(3);
%------iterate to b converge-----------------
p=2;
while   resnorm>1.17e11
    b=b_test(p)
    W0=W0_test(p)
    lamplus=lamplus_test(p)
    lamminus=lamminus_test(p)
    for  i=1:length(NF)
        n=1;
        tensor = [stressben(i)*sind(n*360/stepnumber)+m(i) 0 0 ;...
            0 0 0 ;...
            0 0 0 ];
        sigm=m(i);
        scentre=[2*sigm/3            0                0 ;...
            0             -sigm/3                0 ;...
            0                       0      -sigm/3];
        run('Damiter1.m')
        while n<cycles*stepnumber
            tensor = [stressben(i)*sind(n*360/stepnumber)+m(i), 0, 0 ;...
                0, 0, 0 ;...
                0, 0, 0 ; ];
            hydro=1/3*trace(tensor);
            dev1=tensor-hydro*eye(3)-scentre;
            dev11=dev1(1,1); dev12=dev1(1,2); dev13=dev1(1,3);%to give \dot{dev\Sigma}dt
            dev21=dev1(2,1); dev22=dev1(2,2); dev23=dev1(2,3);
            dev31=dev1(3,1); dev32=dev1(3,2); dev33=dev1(3,3);
            tensor = [stressben(i)*sind((n+1)*360/stepnumber)+m(i), 0, 0 ;...
                0, 0,  0 ;...
                0, 0,  0 ; ];
            run('Damiter2.m')
            n=n+1;
        end
        Smax_bt2d(i)=max(Smax);
        alp_bt2d(i)=mean(alp_ref);
        hydroplus(i)=1/sqrt(2)*(max(hydro)-mean(hydro))+mean(hydro);
        hydrominus(i)=1/sqrt(2)*(min(hydro)-mean(hydro))+mean(hydro);
    end
    parameters0=[W0,b,lamplus,lamminus];  %from pth fit
    [parameters_fit,resnorm,exitflag]=lsqnonlin(fun_analytical,parameters0,lb,ub);
    least_square(p)=resnorm;
    resnorm(end)
    p=p+1;
    W0_test(p)=parameters_fit(1);
    b_test(p)=parameters_fit(2);
    lamplus_test(p)=parameters_fit(3);
    lamminus_test(p)=parameters_fit(3);
    if abs(b_test(p)-b)/b<1e-4 %if new fit stucks
        b_test(p)=b_test(p)*randi([110,200],1,1)/100; %--multiplied by 1.1~2
            if b_test(p)>ub(2)
            b_test(p)=1.01;
        end    
    end
    if abs(W0_test(p)-W0)/W0<1e-5 %if new fit stucks
        W0_test(p)=W0_test(p)*randi([110,900],1,1)/100; %--multiplied by 1.1~9
        if W0_test(p)>ub(1)
            W0_test(p)=lb(1);
        end    
    end
end
W0=W0_test(p);
b=b_test(p);
lamplus=lamplus_test(p);
lamminus=lamminus_test(p);
least_square=resnorm(end)
W0
b
lamplus
lamminus
% lamminus=lamratio.*lamplus;
save('AL6082T6.mat','W0','b','lamplus','lamminus','-append');
% %%%%--------------------------------------------------------------------------
parameters=[W0,b,lamplus,lamminus];
NF_num=fun_analytical(parameters)+NF';


%------------plotting-------------------
figure(1);%----SN---
experiments_bt2d=plot(NF,Smax_bt2d,'ko','MarkerSize',12,'LineWidth', 3);
hold on;
MatlabFit_bt2d=plot(NF_num,Smax_bt2d,'r^','MarkerSize',12,'LineWidth', 3);
xlabel NF;
ylabel Smax;
hLegend=legend([experiments_bt2d,MatlabFit_bt2d],...
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

figure(2);
err_bt = loglog(NF,NF_num,'o','MarkerSize',12,'LineWidth', 3,'MarkerEdgeColor',[	139 69 19]/255, 'MarkerFaceColor','none');
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
hLegend=legend([err_bt],...
    'Bending-torsion test on 30AL6082T6(R=-1)',...
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
saveas(gcf,'F:\Git\Anew\figures\AL6082T6_bt2d_sn.png');
figure(2);
saveas(gcf,'F:\Git\Anew\figures\AL6082T6_bt2d_err.png');

sp=actxserver('SAPI.SpVoice');
sp.Speak('done');

