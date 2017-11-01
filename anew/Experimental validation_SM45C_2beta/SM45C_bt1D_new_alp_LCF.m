% steel SM45C data from Jabbado thesis
%-------------------------------------------------bending--------------------------------
clear Smax_tor Smax_ben alp_tor alp_ben NFben_num NFtor_num
load('gaussian.mat');
%vpa(x,289);
%vpa(weight,289);

NFben=[17520.33619
33991.10739
52427.26419
91077.18396
156881.5804
222261.1047
446114.6594
822487.2063
1279413.916
1453321.25
2440359.869
3428114.749
6880791.155
6213809.484
9342857.067
7240667.386]' ;

stressben=1e6.*[632.13256
590.05764
552.01729
529.5389
506.48415
489.76945
466.7147
463.83285
459.2219
463.83285
454.03458
455.18732
450
437.31988
441.93084
424.0634]';%to get Smaxben

%-------------------------------------------------torsion--------------------------------
stresstor=1e6.*[404.11255
394.87734
375.25253
363.13131
354.4733
345.8153
338.31169
331.38528
329.07648
322.15007
];%to get Smaxtor
NFtor=[27957.25096
47749.27586
76193.65712
100000
162305.0764
182806.9151
296704.9032
575635.5425
822487.2063
2203806.358
];

%%------LCF regime------------
indices =find(NFben<7e5);
NFben=NFben(indices);
stressben=stressben(indices);
indices = find(NFtor<7e5);
NFtor=NFtor(indices);
stresstor=stresstor(indices);
%% 

m=zeros(size(NFben));
sigm=m(1);
scentre=[2*sigm/3            0                0 ;...
                 0             -sigm/3                0 ;...
                 0                       0      -sigm/3]; 
for  i=1:length(NFben)
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
        dev1=tensor-hydro(n)*eye(3)-scentre;
        dev11=dev1(1,1); dev12=dev1(1,2); dev13=dev1(1,3);
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
    hydroplus(i)=1/sqrt(2)*(max(hydro)-mean(hydro))+mean(hydro);
    hydrominus(i)=1/sqrt(2)*(min(hydro)-mean(hydro))+mean(hydro);
end


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

