% steel AL6082T6 data from Jabbado thesis
%-------------------------------------------------bending--------------------------------
clear;
load('AL6082T6.mat');  %update 'W0','b'
load('gaussian.mat');
gam=100;
NFben=[160000
] ;
stressben=1e6.*[590 
];
m=0;
sigm=m(1);
for  i=1:length(NFben)
    n=1;      
    tensor = [stressben(i)*sind(n*360/stepnumber)+m(i) 0 0 ;...
        0 0 0 ;...
        0 0 0 ];
    run('Damiter1.m')
    while n<5*stepnumber
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
        figure(1);
        hold on;
         plot(n,Smax(n)/30.*sign(hydro(n)),'ro');
         plot(n,bs(1024),'ks');
         plot(n,bs(600),'m^');
         plot(n,bs(300),'gv');
        n=n+1;
    end
    Smax_ben(i)=max(Smax); %max sqrt of J2,a
    alp_ben(i)=mean(alp_ref);
    hydroplus(i)=1/sqrt(2)*(max(hydro)-mean(hydro))+mean(hydro);
    hydrominus(i)=1/sqrt(2)*(min(hydro)-mean(hydro))+mean(hydro);
end
%after several cycles the backstress accomondates at a positive value,
%then strats plastic shakedown
bs(1024)
bs(600)
bs(300)
