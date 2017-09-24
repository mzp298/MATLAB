% steel 10HNAP data from Jabbado thesis
%-------------------------------------------------bending--------------------------------
load('HNAP.mat');
load('gaussian.mat');
m=0;
NFben=[1.00E+05
2.00E+05
3.00E+05
4.00E+05
5.00E+05
6.00E+05
7.00E+05
8.00E+05
9.00E+05
1.00E+06] ;
stressben=1e6.*[326.69
304.42
292.11
283.68
277.30
272.20
267.96
264.34
261.19
258.41];%to get Smaxben

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

%% 
%-------------------------------------------------torsion--------------------------------
NFtor=NFben;
stresstor=1e6.*[260.70
242.70
232.17
224.70
218.90
214.16
210.16
206.69
203.63
200.90
];%to get Smaxtor
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
    alp_tor(i)=mean(alp);
end
