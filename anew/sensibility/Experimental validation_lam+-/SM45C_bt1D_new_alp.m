% steel SM45C data from Jabbado thesis
%-------------------------------------------------bending--------------------------------
load('SM45C.mat');
load('gaussian.mat');
%vpa(x,289);
%vpa(weight,289);

m=0;%mean tension
% NFben=[1.9e4 3.6e4 5.8e4 9e4 1.7E5 2.2E5 4.4E5 8E5 1.4E6 1.5E6 2.4E6 3.3E6]' ;
% stressben=[6.2e8 5.9e8 5.52e8 5.35e8 5.05E8 4.9E8 4.70E8 4.65E8 4.62E8 4.65E8 4.6E8 4.58E8 ]';%to get Smaxben
NFben=[1.7E5 2.2E5 4.4E5 8E5 1.4E6 1.5E6 2.4E6 3.3E6]' ;
stressben=[5.05E8 4.9E8 4.70E8 4.65E8 4.62E8 4.65E8 4.6E8 4.58E8 ]';%to get Smaxben

for  i=1:length(NFben)
    n=1;       %initial recording point
    tensor = [stressben(i)*sind(n*360/stepnumber)+m 0 0 ;...
        0 0 0 ;...
        0 0 0 ];
    run('Damiter1.m')
    while n<cycles*stepnumber
        tensor = [stressben(i)*sind(n*360/stepnumber)+m, 0, 0 ;...
            0, 0, 0 ;...
            0, 0, 0 ; ];
        hydro(n)=1/3*trace(tensor);
        dev1=tensor-hydro(n)*eye(3);
        dev11=dev1(1,1); dev12=dev1(1,2); dev13=dev1(1,3);
        dev21=dev1(2,1); dev22=dev1(2,2); dev23=dev1(2,3);
        dev31=dev1(3,1); dev32=dev1(3,2); dev33=dev1(3,3);
        tensor = [stressben(i)*sind((n+1)*360/stepnumber)+m, 0, 0 ;...
            0, 0,  0 ;...
            0, 0,  0 ; ];
        run('Damiter2.m')
         n=n+1;
    end
    Smax_ben(i)=(max(Smax)-min(Smax))/2; %sqrt(2/3).*stressben
    %     Wcyc1(i)=4*(E-k)*(1+nu)*(b-1)/(E*(E+k*nu)*b*(b+1)).*Smax_ben(i).^(b+1).*(y-1/3*lam*m(i)).^(1-b) ;
    alp_ben(i)=mean(alp);
    hydroplus(i)=mean(hydro(find(hydro>0)));
    hydrominus(i)=mean(hydro(find(hydro<0)));
end

%-------------------------------------------------torsion--------------------------------
% stresstor=[4.05E8 
% 3.99E8 
% 3.80E8 
% 3.65E8 
% 3.55E8 
% 3.49E8 
% 3.40E8 
% 3.35E8 
% 3.33E8 
% 3.25E8];%to get Smaxtor
% NFtor=[2.9E4 
% 5E4 
% 7.8E4 
% 1E5 
% 1.8E5 
% 1.9E5 
% 3E5 
% 5.8E5 
% 8.2E5 
% 2.25E6];
stresstor=[
3.65E8 
3.55E8 
3.49E8 
3.40E8 
3.35E8 
3.33E8 
3.25E8];%to get Smaxtor
NFtor=[
1E5 
1.8E5 
1.9E5 
3E5 
5.8E5 
8.2E5 
2.25E6];

for  i=1:length(NFtor)
    n=1;       %initial recording point
    tensor = [0 stresstor(i)*sind(n*360/stepnumber) 0 ;...
        stresstor(i)*sind(n*360/stepnumber) 0  0 ;...
        0 0 0 ];
     run('Damiter1.m')
      
    while n<cycles*stepnumber
        tensor = [0, stresstor(i)*sind(n*360/stepnumber), 0 ;...
            stresstor(i)*sind(n*360/stepnumber), 0,  0 ;...
            0,0, 0; ];
        hydro=1/3*trace(tensor);
        dev1=tensor-hydro*eye(3);
        dev11=dev1(1,1); dev12=dev1(1,2); dev13=dev1(1,3);%to give \dot{dev\Sigma}dt
        dev21=dev1(2,1); dev22=dev1(2,2); dev23=dev1(2,3);
        dev31=dev1(3,1); dev32=dev1(3,2); dev33=dev1(3,3);
        tensor = [0, stresstor(i)*sind((n+1)*360/stepnumber), 0 ;...
            stresstor(i)*sind((n+1)*360/stepnumber), 0,  0 ;...
            0, 0, 0; ];
                run('Damiter2.m')
        n=n+1;
    end
    Smax_tor(i)=(max(Smax)-min(Smax))/2; %sqrt(2/3)*stresstor
    %     Wcyc1(i)=4*(E-k)*(1+nu)*(b-1)/(E*(E+k*nu)*b*(b+1)).*Smax_tor(i).^(b+1).*(y-1/3*lam*m(i)).^(1-b) ;
    alp_tor(i)=mean(alp);
end

% figure(5)
% semilogx(NFben,stressben,'*b')
% hold on
% semilogx(NFtor,stresstor,'or')
% m=196e6;%mean tension
% NFben_m=[4.50E+04
% 5.40E+04
% 7.00E+04
% 7.10E+04
% 1.10E+05
% 1.50E+05
% 2.00E+05
% 2.10E+05
% 3.00E+05
% 4.30E+05
% 6.90E+05
% 5.80E+05] ;
% stressben_m=1e6.*[540
% 515
% 520
% 485
% 485
% 475
% 460
% 455
% 435
% 415
% 410
% 390];%to get Smaxben
% semilogx(NFben_m,stressben_m,'sk')