function [newX,idkeep] = transform2isotrop(oldX, factor)
%==============================================
% to transform any array into an isotrope array
% synopsis
%     newX = transform2isotrop(oldX)
% oldX : array M x D
%     where M is the number of sensors
%           D the space dimension, ie 2 or 3
% factor : radius factor w.r.t. the highest
%          eigenvalue
%==============================================
%
[M,D]         = size(oldX);
distold       = sum(oldX .^2,2);
[ss,idkeep]   = sort(distold);
oldX          = oldX(idkeep(1:M-D),:);
XXT           = oldX'*oldX;
[UU,DD]       = eig(XXT);
DD            = diag(DD);
% any number greater than max(DD)
% provides a solution
MMM           = max(DD)*factor;
xMp           = zeros(D,D);
for id = 1:D
    xMp(:,id) = sqrt(MMM-DD(id))*UU(:,id);
end
newX          = [oldX ; xMp'];
%==============================================
