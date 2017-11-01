%-------------repeated code to get integration on all scales----------
hydro(n+1)=1/3*trace(tensor);
if hydro(n+1)>=0
yield(n+1)=y-lamplus(end)*hydro(n+1); %yield stress at time step n+1
else
yield(n+1)=y-lamminus(end)*hydro(n+1);    
end
scentre=[2*sigm/3            0                0 ;...
                 0             -sigm/3                0 ;...
                 0                       0      -sigm/3];
devn=tensor-hydro(n+1)*eye(3)-scentre;
dev11g=devn(1,1); dev12g=devn(1,2); dev13g=devn(1,3);
dev21g=devn(2,1); dev22g=devn(2,2); dev23g=devn(2,3);
dev31g=devn(3,1); dev32g=devn(3,2); dev33g=devn(3,3);
Smax(n+1)=1/sqrt(2).*norm(devn,'fro');


trial11=bsxfun(@plus,Sb11,(dev11g-dev11)); trial12=bsxfun(@plus,Sb12,(dev12g-dev12));trial13=bsxfun(@plus,Sb13,(dev13g-dev13));
trial21=bsxfun(@plus,Sb21,(dev21g-dev21)); trial22=bsxfun(@plus,Sb22,(dev22g-dev22));trial23=bsxfun(@plus,Sb23,(dev23g-dev23));
trial31=bsxfun(@plus,Sb31,(dev31g-dev31)); trial32=bsxfun(@plus,Sb32,(dev32g-dev32));trial33=bsxfun(@plus,Sb33,(dev33g-dev33));

trialtensor=[trial11; trial12; trial13; trial21; trial22; trial23;trial31; trial32; trial33];
% Smaxtrial=1/sqrt(2).*sqrt(sum(trialtensor.^2));
Smaxtrial=sqrt(trial11.^2+trial12.^2);
eta=bsxfun(@minus,bsxfun(@times,Smaxtrial/yield(n+1),s),1); %1*64
eta(eta<0)=0;

Sb11=bsxfun(@rdivide,trial11,bsxfun(@plus,eta,1));Sb12=bsxfun(@rdivide,trial12,bsxfun(@plus,eta,1));Sb13=bsxfun(@rdivide,trial13,bsxfun(@plus,eta,1));
Sb21=bsxfun(@rdivide,trial21,bsxfun(@plus,eta,1));Sb22=bsxfun(@rdivide,trial22,bsxfun(@plus,eta,1));Sb23=bsxfun(@rdivide,trial23,bsxfun(@plus,eta,1));
Sb31=bsxfun(@rdivide,trial31,bsxfun(@plus,eta,1));Sb32=bsxfun(@rdivide,trial32,bsxfun(@plus,eta,1));Sb33=bsxfun(@rdivide,trial33,bsxfun(@plus,eta,1));

Ws=(bsxfun(@minus,Smaxtrial,bsxfun(@rdivide, yield(n+1),s))<=0).*...
    (0)+...
    (bsxfun(@minus,Smaxtrial,bsxfun(@rdivide, yield(n+1),s))>0).*...
    ((E-k)*(1+nu)*(2*E*(E+k*nu))^-1*bsxfun(@times,weight,bsxfun(@rdivide,bsxfun(@times,bsxfun(@minus,Smaxtrial,bsxfun(@rdivide, yield(n+1),s)),yield(n+1)),s)));
W= sum(Ws);
W_accumulate(n+1)=W_accumulate(n)+W;
sequence=((Smax(n+1)*yield(n+1)^-1)*(1-Smax(n+1)*yield(n+1)^-1)^-1)^fb;
sequence(sequence<0)=0;
alp(n+1)=1-a*sequence;


if n+1>(cycles-1)*stepnumber% last cycle after adaptation
alp_ref(1)=alp(n);
n_ref(1)=n;
W_ref(1)=W;
if abs(alp(n+1)-alp_ref(g))>delta_alp %----giving scalar value to iteration after the addaptation cycle(decrease time step)
        alp_ref(g+1)=alp(n+1);
        n_ref(g+1)=n+1;
        W_ref(g+1)=W_accumulate(n+1)-W_accumulate(n_ref(g));
    g=g+1;
end
end


%     if abs(alp(n+1)-alp_ref(g))>delta_alp %----giving scalar value to iteration after the addaptation cycle(decrease time step)
%         alp_ref(g+1)=alp(n);
%         n_ref(g+1)=n;
%         W_ref(g+1)=W_accumulate(n+1)-W_accumulate(n_ref(g));
%     g=g+1;
%     end
