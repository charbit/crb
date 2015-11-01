%===================================================
function [area] = ellipse(X0, R, alpha)
%===================================================
% Drawing ellipse
% SYNOPSIS: ELLIPSE(X0, R, alpha)
% inputs:
%    X0 = Coordinates of the ellipse's center (2x1)
%    R  = A covariance positive (2x2) matrix
%    alpha  = confidence level 
% outputs:
% areaMC : area using monte-carlo
% area : area of the confidence ellipse at 100alpha%
%===================================================
c=-2*log(1-alpha);
N=100; theta = (0:N) * (2*pi) ./ N ;
Y = sqrt(c)*[cos(theta);sin(theta)];
Fm1=sqrtm(R);
X = diag(X0)*ones(2,N+1)+Fm1*Y;
plot(X(1,:),X(2,:));
grid on
valp   = eig(R);
area   = sqrt(prod(valp))*c*pi;
% %=====
% G=100000;
% invR=inv(R);
% AA=sqrt(max(valp*c));
% points=2*AA*(rand(G,2)-1/2);
% ss=zeros(G,1);
% for ii=1:G
%     ppaux=points(ii,:);
%     ss(ii)=ppaux*invR*ppaux'<c;
% end
% areaMC = sum(ss)*4*AA*AA/G;
%=====================================================