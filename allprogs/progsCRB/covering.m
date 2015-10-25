clear
M=10;
Lruns=1000;
mr = zeros(Lruns,1);
x=zeros(M,2,Lruns);
d=zeros(M*(M-1)/2,Lruns);
UU=1;
for ir=1:Lruns
    if UU
        rho = rand(M,1);
        phi = 360*rand(M,1);
    else
        rho = ones(M,1);
        phi = (0:M-1)'*360/M;
    end
    x(:,1,ir) = rho .* cos(phi*pi/180);
    x(:,2,ir) = rho .* sin(phi*pi/180);
    cp=0;
    for im1=1:M-1
        for im2=im1+1:M
            cp=cp+1;
            d(cp,ir)=norm(x(im1,:,ir)-x(im2,:,ir));
        end
    end
    mr(ir) = min(d(:,ir));
end

[mm,indmm] = max(mr);
%%
figure(1)
if UU
    subplot(223)
    
    plot(x(:,1,indmm),x(:,2,indmm),'o')
    grid on
    set(gca,'xlim',[-1 1])
    set(gca,'ylim',[-1 1])
    axis('square')
    title(sprintf('dmin = %5.2f',mm))
    
    subplot(224)
    plot(sort(d(:,indmm)),'.')
    grid on
    
else
    subplot(221)
    
    plot(x(:,1,indmm),x(:,2,indmm),'o')
    grid on
    set(gca,'xlim',[-1 1])
    set(gca,'ylim',[-1 1])
    axis('square')
    title(sprintf('dmin = %5.2f',mm))
    
    subplot(222)
    plot(sort(d(:,indmm)),'.')
    grid on
end
%%

newX = anyarray2isotrop(squeeze(x(:,:,indmm))');

XXTnew = newX*newX'

%%
HorizontalSize = 16;
VerticalSize   = 10;
set(gcf,'units','centimeters');
set(gcf,'paperunits','centimeters');
set(gcf,'PaperType','a3');
set(gcf,'position',[0 5 HorizontalSize VerticalSize]);
%         set(gcf,'position',[0 5 HorizontalSize VerticalSize]);
set(gcf,'paperposition',[0 0 HorizontalSize VerticalSize]);
set(gcf,'color', [1,1,0.92]);
set(gcf, 'InvertHardCopy', 'off');


printdirectory  = ' ../../figures/';
fileprintepscmd = sprintf('print -depsc -loose %suniforinterdistances.eps',printdirectory);
fileeps2pdfcmd  = sprintf('!/Library/TeX/texbin/epstopdf %suniforinterdistances.eps',printdirectory);
filermcmd       = sprintf('!rm %suniforinterdistances.eps',printdirectory);

figure(1)
% eval(fileprintepscmd)
% eval(fileeps2pdfcmd)
% eval(filermcmd)