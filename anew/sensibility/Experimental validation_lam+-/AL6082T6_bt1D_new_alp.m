% steel AL6082T6 data from Jabbado thesis
%-------------------------------------------------bending--------------------------------
load('AL6082T6.mat');  %update 'W0','b'
load('gaussian.mat');
NFben=[160000
248518
444411
1069220 
56285 
1238325 
200480 
423590] ;
stressben=1e6.*[190 
180 
164 
144 
224 
145 
187 
161];
m=zeros(size(NFben));

for  i=1:length(stressben) %experimental points index
    n=1;
    tensor = [stressben(i)*sind(n*360/stepnumber)+m(i) 0 0 ;...
        0 0 0 ;...
        0 0 0 ];
    run('Damiter1.m')
    while n<cycles*stepnumber
        tensor = [stressben(i)*sind(n*360/stepnumber)+m(i), 0, 0 ;...
            0, 0, 0 ;...
            0, 0, 0 ; ];
        hydro(n)=1/3*trace(tensor);
        dev1=tensor-hydro(n)*eye(3);
        dev11=dev1(1,1); dev12=dev1(1,2); dev13=dev1(1,3);
        dev21=dev1(2,1); dev22=dev1(2,2); dev23=dev1(2,3);
        dev31=dev1(3,1); dev32=dev1(3,2); dev33=dev1(3,3);
        tensor = [stressben(i)*sind((n+1)*360/stepnumber)+m(i), 0, 0 ;...
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
NFtor=[534032
    26987
    76665
    132295
    203535
    16195
    1.1E6
    565150
    ];
stresstor=1e6.*[117
    155
    127
    127
    117
    155
    106
    104
    ];%to get Smaxtor
% 
% %------------reduce to HCF regime-----------
% NFtor=[534032
%     76665
%     132295
%     203535
%     1.1E6
%     565150
%     ];
% stresstor=1e6.*[117
%     127
%     127
%     117
%     106
%     104
%     ];%to get Smaxtor

m=zeros(size(NFtor));

for  i=1:length(stresstor) %experimental points index
    n=1;
    tensor = [0 stresstor(i)*sind(n*360/stepnumber)+m(i) 0 ;...
        stresstor(i)*sind(n*360/stepnumber)+m(i) 0 0 ;...
        0 0 0 ];
    run('Damiter1.m')
    while n<cycles*stepnumber
        tensor = [0 stresstor(i)*sind(n*360/stepnumber)+m(i) 0 ;...
            stresstor(i)*sind(n*360/stepnumber)+m(i) 0 0 ;...
            0  0  0 ; ];
        hydro=1/3*trace(tensor);
        dev1=tensor-hydro*eye(3);
        dev11=dev1(1,1); dev12=dev1(1,2); dev13=dev1(1,3);
        dev21=dev1(2,1); dev22=dev1(2,2); dev23=dev1(2,3);
        dev31=dev1(3,1); dev32=dev1(3,2); dev33=dev1(3,3);
        tensor = [0 stresstor(i)*sind((n+1)*360/stepnumber)+m(i) 0 ;...
            stresstor(i)*sind((n+1)*360/stepnumber)+m(i) 0   0 ;...
            0  0   0 ; ];
        run('Damiter2.m')
        n=n+1;
    end
    Smax_tor(i)=(max(Smax)-min(Smax))/2; %sqrt(2/3)*stresstor
    alp_tor(i)=mean(alp);
end
