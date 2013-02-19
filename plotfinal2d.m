function [hLa] = plotfinal2d(name,nc,res,xm,ym,subfiglabel)

x = nc_varget(nc,['x' res]); 
y = nc_varget(nc,['y' res]); 
tmp = nc_varget(nc,['Q' res]);
q0 = squeeze(tmp(1,:,:));
qf = squeeze(tmp(end,:,:));

cntr = [-0.05 -0.01 0.05:0.1:0.95];

hLa = axes('Position',[xm ym 0.213 0.288]);

      cntr2 = [-0.55:0.1:-0.05 -0.01 0.05:0.3:0.95 1.01 1.05:0.1:1.55];
      contourf(x,y,qf,cntr2); caxis([-0.5 1.5]);
      colormap(blankstep); hold on;
% $$$       [c,h] = v6contour(x,y,q0,[0.5 0.5],'k'); 
% $$$       set(h,'LineWidth', 1.25);
      axis square
      set(hLa,'XTickLabel','');
      set(hLa,'YTickLabel','');

clear hLt
tmp = min(min(qf));
if tmp<0 & tmp>-1e-3
  hLt(1) = text(0.55,0.17,sprintf('Min= %6.1e',min(min(qf))));
else
  hLt(1) = text(0.55,0.17,sprintf('Min= %5.3f',min(min(qf))));
end
hLt(2) = text(0.55,0.09,sprintf('Max= %5.3f',max(max(qf))));
hLt(3) = text(0.05,0.17,['E_2= ' sprintf('%5.3f', ...
                                         sqrt(mean((qf(:)-q0(:)).^2)))]);
hLt(4) = text(0.05,0.09,['E_{\infty}= ' sprintf('%5.3f', ...
                                                max(abs(qf(:)-q0(:))))]);
set(hLt,'VerticalAlignment','Cap','FontSize',8);

hLu = text(0.02,0.93,[subfiglabel ') ' name]); set(hLu,'FontSize', ...
                                                       8);
set(hLu,'VerticalAlignment','Cap');
set(hLa,'XTickLabel','');
set(hLa,'YTickLabel','');

