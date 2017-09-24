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
m=0;
sigm=m;
scentre=[2*sigm/3            0                0 ;...
                 0             -sigm/3                0 ;...
                 0                       0      -sigm/3]; 
for  i=1:length(NFben)
    n=1;      
    tensor = [stressben(i)*sind(n*360/stepnumber)+m 0 0 ;...
        0 0 0 ;...
        0 0 0 ];
    run('Damiter1.m')
    while n<cycles*stepnumber
        tensor = [stressben(i)*sind(n*360/stepnumber)+m, 0, 0 ;...
            0, 0, 0 ;...
            0, 0, 0 ; ];
        hydro(n)=1/3*trace(tensor);
        dev1=tensor-hydro(n)*eye(3)-scentre;
        dev11=dev1(1,1); dev12=dev1(1,2); dev13=dev1(1,3);
        dev21=dev1(2,1); dev22=dev1(2,2); dev23=dev1(2,3);
        dev31=dev1(3,1); dev32=dev1(3,2); dev33=dev1(3,3);
        tensor = [stressben(i)*sind((n+1)*360/stepnumber)+m, 0, 0 ;...
            0, 0,  0 ;...
            0, 0,  0 ; ];
        run('Damiter2.m')
        n=n+1;
    end
    Smax_ben(i)=max(Smax); %max sqrt of J2,a
    alp_ben(i)=mean(alp);
    hydroplus(i)=1/sqrt(2)*(max(hydro)-mean(hydro))+mean(hydro);
    hydrominus(i)=1/sqrt(2)*(min(hydro)-mean(hydro))+mean(hydro);
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
    alp_tor(i)=mean(alp);
end
