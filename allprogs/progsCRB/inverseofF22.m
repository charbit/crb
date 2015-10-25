% clear
K=12;
% M33 = randn(3);
% F123 = M33*M33';
% ZK = zeros(3,K);
% Zsigma=zeros(3,1);
% fsigma2=0.2;
unkn = randn(1,K);

fsigma2=20;
MKK =randn(K);
FKK = MKK*MKK';

F = [F123 Zsigma ZK ; ...
    Zsigma' 1/fsigma2 unkn;...
    ZK' unkn' FKK];

Finv = inv(F);

CRB1 = Finv(1:3,1:3);

[CRB1 inv(F123)]