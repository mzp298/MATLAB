% %%------------------1 torsion fit only--------------
fun_analytical=@(parameters)parameters(1).*E.*(E+k.*nu).*parameters(2).*(parameters(2)+1).*...
    (2.*(1-[alp_tor]).*(E-k).*(1+nu).*(parameters(2)-1)).^-1.*...
    ([Smax_tor].^(parameters(2)+1).*(y-lamplus.*[zeros(size(stresstor))]').^(1-parameters(2))...
    +[Smax_tor].^(parameters(2)+1).*(y-lamminus.*[zeros(size(stresstor))]').^(1-parameters(2))).^-1-[NFtor]';
% %-------------lsqnonlin locally stuck, issue new parameters of the non-continuous function using random
% %-------------generated new values-------------
parameters=[3e8,6];
lb  =  [0,   1];
ub =  [1e12, ub_beta];    %bmax=1 gives 7.26e+06
[parameters_fit,resnorm,exitflag]=lsqnonlin(fun_analytical,parameters,lb,ub);
[W0,b]=deal(parameters_fit(1),parameters_fit(2));
least_square=resnorm;
W0_test(2)=parameters_fit(1); %original fit
b_test(2)=parameters_fit(2); %original fit

%------iterate to b converge-----------------
p=2;
while   resnorm>resnorm_limit
    b=b_test(p)
    W0=W0_test(p)
 for  i=1:length(stresstor) %experimental points index
    n=1;
    tensor = [0 stresstor(i)*sind(n*360/stepnumber) 0 ;...
        stresstor(i)*sind(n*360/stepnumber) 0  0 ;...
        0 0 0 ];
    run('Damiter1.m')
    while n<cycles*stepnumber
        tensor = [0, stresstor(i)*sind(n*360/stepnumber), 0 ;...
            stresstor(i)*sind(n*360/stepnumber), 0,  0 ;...
            0,0, 0; ];
        hydro=1/3*trace(tensor);
        dev1=tensor-hydro*eye(3)-scentre; %-----------scentre should be 0 in torsion
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
    alp_tor(i)=mean(alp_ref);
end
    parameters0=[W0,b];  %from pth fit
    [parameters_fit,resnorm,exitflag]=lsqnonlin(fun_analytical,parameters0,lb,ub);
    least_square(p)=resnorm;
    resnorm(end)
    p=p+1;
    W0_test(p)=parameters_fit(1);
    b_test(p)=parameters_fit(2);
    if abs(b_test(p)-b)/b<1e-2 %if new fit stucks
        b_test(p)=b_test(p)*randi([110,200],1,1)/100; %--multiplied by 1.1~2
        if b_test(p)>ub(2)
            b_test(p)=1.01;
        end
    end
    if abs(W0_test(p)-W0)/W0<1e-2 %if new fit stucks
        W0_test(p)=W0_test(p)*randi([110,900],1,1)/100; %--multiplied by 1.1~9
        if W0_test(p)>ub(1)
            W0_test(p)=lb(1);
        end
    end
end
W0=W0_test(p);
b=b_test(p);

least_square=resnorm(end)
W0
b
lamplus
lamminus

parameters=[W0,b];
NFtor_num=fun_analytical(parameters)+[NFtor]';

%%

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
    Smax_ben(i)=max(Smax); %max sqrt of J2,a
    alp_ben(i)=mean(alp_ref);
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
% %%
% 
%%----------3 num_opt to get NFtor_num------------------------
    for  i=1:length(stresstor) %experimental points index
    n=1;
    tensor = [0 stresstor(i)*sind(n*360/stepnumber) 0 ;...
        stresstor(i)*sind(n*360/stepnumber) 0  0 ;...
        0 0 0 ];
    run('Damiter1.m')
    while n<cycles*stepnumber
        tensor = [0, stresstor(i)*sind(n*360/stepnumber), 0 ;...
            stresstor(i)*sind(n*360/stepnumber), 0,  0 ;...
            0,0, 0; ];
        hydro=1/3*trace(tensor);
        dev1=tensor-hydro*eye(3)-scentre; %-----------scentre should be 0 in torsion
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
    alp_tor(i)=mean(alp_ref);
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
    NFtor_num(i)=e/stepnumber
    end
save('SM45C.mat','NFben','NFtor','NFben_num','NFtor_num','lamplus','lamminus','-append');

% plot(Smax,'rs')
% hold on
% plot(yield/max(s),'ko')

