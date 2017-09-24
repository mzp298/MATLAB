%%------------------1 torsion fit only--------------
fun_analytical=@(parameters)parameters(1).*E.*(E+k.*nu).*parameters(2).*(parameters(2)+1).*...
    (2.*(1-[alp_tor]).*(E-k).*(1+nu).*(parameters(2)-1)).^-1.*...
    ([Smax_tor].^(parameters(2)+1).*(y-parameters(3).*[zeros(size(stresstor))]').^(1-parameters(2))...
    +[Smax_tor].^(parameters(2)+1).*(y-parameters(4).*[zeros(size(stresstor))]').^(1-parameters(2))).^-1-[0.65.*NFtor]';
% %-------------lsqnonlin locally stuck, issue new parameters of the non-continuous function using random
% %-------------generated new values-------------
parameters=[1e5,3,0,0];
lb  =  [1e3,   1,0,0];
ub =  [1e12,   100,5,5];
[parameters_fit,resnorm,exitflag]=lsqnonlin(fun_analytical,parameters,lb,ub);
[W0,b,lamplus,lamminus]=deal(parameters_fit(1),parameters_fit(2),parameters_fit(3),parameters_fit(4));
least_square=resnorm;
W0_test(2)=parameters_fit(1); %original fit
b_test(2)=parameters_fit(2); %original fit
lamplus_test(2)=parameters_fit(3);
lamminus_test(2)=parameters_fit(4);
%------iterate to b converge-----------------
p=2;
while   resnorm>2.1e12
    b=b_test(p)
    W0=W0_test(p)
    lamplus=lamplus_test(p)
    lamminus=lamminus_test(p)
    for  i=1:length(NFben)
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
            hydro(n)=1/3*trace(tensor);
            dev1=tensor-hydro(n)*eye(3)-scentre;
            dev11=dev1(1,1); dev12=dev1(1,2); dev13=dev1(1,3);%to give \dot{dev\Sigma}dt
            dev21=dev1(2,1); dev22=dev1(2,2); dev23=dev1(2,3);
            dev31=dev1(3,1); dev32=dev1(3,2); dev33=dev1(3,3);
            tensor = [stressben(i)*sind((n+1)*360/stepnumber)+m(i), 0, 0 ;...
                0, 0,  0 ;...
                0, 0,  0 ; ];
            run('Damiter2.m')
            n=n+1;
        end
        Smax_ben(i)=max(Smax);
        alp_ben(i)=mean(alp);
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
    lamminus_test(p)=parameters_fit(4);
    if abs(b_test(p)-b)/b<1e-5 %if new fit stucks
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
save('SM45C.mat','W0','b','lamplus','lamminus','-append');
parameters=[W0,b,lamplus,lamminus];
NFtor_num=fun_analytical(parameters)+[0.65.*NFtor]';


%%-------------2 analytical of NF_ben to identify lam+-
fun_analytical_ben=@(parameters)W0.*E.*(E+k.*nu).*b.*(b+1).*...
    (2.*(1-[alp_ben]).*(E-k).*(1+nu).*(b-1)).^-1.*...
    ([Smax_ben].^(b+1).*(y-parameters(1).*[hydroplus]).^(1-b)...
    +[Smax_ben].^(b+1).*(y-parameters(2).*[hydrominus]).^(1-b)).^-1-[NFben]';
% parameters=[0.1,0.1];
% lb  =  [0,0];
% ub =  [20,20];
% [parameters_fit,resnorm,exitflag]=lsqnonlin(fun_analytical_ben,parameters,lb,ub);
% [lamplus,lamminus]=deal(parameters_fit(1),parameters_fit(2));
% least_square=resnorm;
% lamplus_test(2)=parameters_fit(1);
% lamminus_test(2)=parameters_fit(2);
% %------iterate to b converge-----------------
% p=2;
% while  abs(lamplus_test(p)-lamplus_test(p-1))/lamplus>1e-2
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
%             hydro(n)=1/3*trace(tensor);
%             dev1=tensor-hydro(n)*eye(3)-scentre;
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
%         alp_ben(i)=mean(alp);
%         hydroplus(i)=1/sqrt(2)*(max(hydro)-mean(hydro))+mean(hydro);
%         hydrominus(i)=1/sqrt(2)*(min(hydro)-mean(hydro))+mean(hydro);
%     end
%     parameters0=[lamplus,lamminus];  %from pth fit
%     [parameters_fit,resnorm,exitflag]=lsqnonlin(fun_analytical_ben,parameters0,lb,ub);
%     least_square(p)=resnorm;
%     resnorm(end)
%     p=p+1;
%     lamplus_test(p)=parameters_fit(1)
%     lamminus_test(p)=parameters_fit(2);
% end
% lamplus=lamplus_test(p);
% lamminus=lamminus_test(p);
% least_square=resnorm(end)
% lamplus
% lamminus
% save('SM45C.mat','lamplus','lamminus','-append');
parameters=[lamplus,lamminus];
NFben_ana=fun_analytical_ben(parameters)+[NFben]';

lamplus=lamplus_num %numerical method is more influenced by lamplus
%%----------3 num_opt to get NFben_num------------------------
for  i=1:length(NFben)
    n=1;
    tensor = [stressben(i)*sind(n*360/stepnumber)+m(i) 0 0 ;...
        0 0 0 ;...
        0 0 0 ];
    sigm=m(i);
    run('Damiter1.m')
    g=1; %first reference index when creating
    while n<cycles*stepnumber
        tensor = [stressben(i)*sind(n*360/stepnumber)+m(i), 0, 0 ;...
            0, 0, 0 ;...
            0, 0, 0 ; ];
        hydro(n)=1/3*trace(tensor);
        dev1=tensor-hydro(n)*eye(3)-scentre;
        dev11=dev1(1,1); dev12=dev1(1,2); dev13=dev1(1,3);%to give \dot{dev\Sigma}dt
        dev21=dev1(2,1); dev22=dev1(2,2); dev23=dev1(2,3);
        dev31=dev1(3,1); dev32=dev1(3,2); dev33=dev1(3,3);
        tensor = [stressben(i)*sind((n+1)*360/stepnumber)+m(i), 0, 0 ;...
            0, 0,  0 ;...
            0, 0,  0 ; ];
        run('Damiter2.m')
        n=n+1;
    end
    Smax_ben(i)=max(Smax);
    alp_ben(i)=mean(alp);
    hydroplus(i)=1/sqrt(2)*(max(hydro)-mean(hydro))+mean(hydro);
    hydrominus(i)=1/sqrt(2)*(min(hydro)-mean(hydro))+mean(hydro);
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
    NFben_num(i)=e/stepnumber
end
save('SM45C.mat','NFben_num','NFben_ana','NFtor_num','alp_ben','alp_tor','lamplus','lamminus','-append');

