function [newX,idkeep] = transform2isotrop(oldX)
%==============================================
% to transform any array into an isotrope array
% synopsis
%     newX = transform2isotrop(oldX)
% oldX : array M x d
%     where M is the number of sensors
%           d the space dimension, ie 2 or 3
%==============================================
%
[M,d]   = size(oldX);
indperm = randperm(M);
distorig = sum(oldX .^2,2);
[ss,idkeep]  = sort(distorig);
oldX    = oldX(idkeep(1:M-d),:);
XXT     = oldX'*oldX;
[UU,DD] = eig(XXT);
DD      = diag(DD);
% any number greater than max(DD)
% provides a solution
MMM     = max(DD)*1.2;
xMp     = zeros(d,d);
for id=1:d
    xMp(:,id) = -sqrt(MMM-DD(id))*UU(:,id);
end
newX = [oldX ; xMp'];
%==============================================
