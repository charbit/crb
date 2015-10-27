function newX = transform2isotrop(gX)
% to transform in isotrope array
[d,M]   = size(gX);
XXT     = gX*gX';
[UU,DD] = eig(XXT);
DD      = diag(DD);
MMM     = max(DD);
xMp = zeros(d,d);
for id=1:d
    xMp(:,id) = sqrt(MMM-DD(id))*UU(:,id);
end
newX = [gX xMp];

