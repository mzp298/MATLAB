%-----------------------2 analytical original(without logform)------------------------
fun_analytical=@(parameters)parameters(1).*E.*(E+k.*nu).*parameters(2).*(parameters(2)+1).*...
    (2.*(1-[alp_tor]).*(E-k).*(1+nu).*(parameters(2)-1)).^-1.*...
    ([Smax_tor].^(parameters(2)+1).*(y-parameters(3).*[zeros(size(stresstor))]').^(1-parameters(2))...
    +[Smax_tor].^(parameters(2)+1).*(y-parameters(3).*[zeros(size(stresstor))]').^(1-parameters(2))).^-1-[NFtor]';

% %%-------------lsqnonlin locally stuck, issue new parameters of the non-continuous function using random
% %%-------------generated new values-------------
parameters=[1e5,1.1,lamplus_initial];
lb  =  [1e3,   1,0];
ub =  [1e10,   10,2];
[parameters_fit,resnorm,exitflag]=lsqnonlin(fun_analytical,parameters,lb,ub);
[W0,b,lamplus,lamminus]=deal(parameters_fit(1),parameters_fit(2),parameters_fit(3),parameters_fit(3));
least_square=resnorm;
W0_test(2)=parameters_fit(1); %original fit
b_test(2)=parameters_fit(2); %original fit
lamplus_test(2)=parameters_fit(3);
lamminus_test(2)=parameters_fit(3);
%------iterate to b converge-----------------
p=2;
while   resnorm>3e8
    b=b_test(p)
    W0=W0_test(p)
    lamplus=lamplus_test(p)
    lamminus=lamminus_test(p)
    for  i=1:length(NFtor)
        n=1;
    tensor = [0 stresstor(i)*sind(n*360/stepnumber) 0 ;...
        stresstor(i)*sind(n*360/stepnumber) 0  0 ;...
        0 0 0 ];
        sigm=m(i);
        scentre=[0            0                0 ;...
            0             0                0 ;...
            0                       0      0];
        run('Damiter1.m')
        while n<cycles*stepnumber
        tensor = [0, stresstor(i)*sind(n*360/stepnumber), 0 ;...
            stresstor(i)*sind(n*360/stepnumber), 0,  0 ;...
            0,0, 0; ];
            hydro=1/3*trace(tensor);
            dev1=tensor-hydro*eye(3)-scentre;
            dev11=dev1(1,1); dev12=dev1(1,2); dev13=dev1(1,3);%to give \dot{dev\Sigma}dt
            dev21=dev1(2,1); dev22=dev1(2,2); dev23=dev1(2,3);
            dev31=dev1(3,1); dev32=dev1(3,2); dev33=dev1(3,3);
        tensor = [0, stresstor(i)*sind((n+1)*360/stepnumber), 0 ;...
            stresstor(i)*sind((n+1)*360/stepnumber), 0,  0 ;...
            0, 0, 0; ];
            run('Damiter2.m')
            n=n+1;
        end
    Smax_tor(i)=max(Smax); %max sqrt of J2,a
    alp_tor(i)=mean(alp);
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
lamminus=lamratio.*lamplus;
save('NCD16.mat','W0','b','lamplus','lamminus','-append');
parameters=[W0,b,lamplus];

NFtor_num=fun_analytical(parameters)+[NFtor]';
NFben_num=parameters(1).*E.*(E+k.*nu).*parameters(2).*(parameters(2)+1).*...
    (2.*(1-[alp_ben]).*(E-k).*(1+nu).*(parameters(2)-1)).^-1.*...
    ([Smax_ben].^(parameters(2)+1).*(y-parameters(3).*[hydroplus]).^(1-parameters(2))...
    +[Smax_ben].^(parameters(2)+1).*(y-parameters(3).*[hydrominus]).^(1-parameters(2))).^-1;

save('NCD16.mat','NFben_num','NFtor_num','alp_ben','alp_tor','-append'); %save the best fit W0, b, alp_ben, alp_tor
