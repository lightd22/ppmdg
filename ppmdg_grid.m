% ----
% Grid testing for PPMDG Advection
% By: Devin Light
% ----

nx = 25;
ny = 25;

% X-grid
% ---
num_elem = nx/5;
elemdx = 1.0/num_elem;
ecent = zeros(1,num_elem);

% Nodal locations in [-1,1]
node4 = [-1.0, -0.654653670707978,0.0, 0.654653670707977,1.0];

ecent(1) = 0.5*elemdx;

 for i = 2:num_elem
   ecent(i)=ecent(i-1)+elemdx;
 end

 x = zeros(1,nx);
 for i = 0:num_elem-1
    x(1+i*5:1+i*5+4) = ecent(i+1)+0.5*elemdx*node4(1:5); 
 end

 xf = zeros(1,nx+1);
 xf(1) = 0.0;
 xf(nx+1) = 1.0;
 for i = 1:nx-1
    xf(i+1) = x(i) + 0.5d0*(x(i+1)-x(i));
 end

 dx = zeros(1,nx);
 for i = 2:nx+1
    dx(i) = xf(i)-xf(i-1);
 end
 
 % Y-grid
 % ---
 
 dy = 1.0/ny;
 
 yf = zeros(1,ny+1);
 y = zeros(1,ny);
 
 yf(1) = 0.0;
 
for i = 1:ny
    yf(i+1) = yf(i) + dy;
    y(i) = 0.5*(yf(i+1)+yf(i));
end

% Even spaced X-grid
% ---
 
xe = zeros(1,nx);
dxe = 1.0/nx;

xe(1) = dxe/2;
for i = 1:nx-1
    xe(i+1) = xe(i)+dxe;
end
xef = zeros(1,nx+1);
xef(1) = 0.0;
for i = 1:nx
    xef(i+1)= xef(i)+ dxe;
end

 [X,Y] = meshgrid(x,y);
 [XF,YF] = meshgrid(xf,yf);
 [XE,YE] = meshgrid(xe,y);
 [ECENT,Y2] = meshgrid(ecent,y);
 [XEF,YEF] = meshgrid(xef,yf);
 hold on
 plot(X,Y,'r.')
 %plot(XF,YF,'rx')
 plot(XE,YE,'ko')
 plot(ECENT,Y2,'g.')
 plot(XEF,YEF,'co')
 for i = 1:ny+1
    line([0 1],[yf(i) yf(i)],'Marker','.','LineStyle','--')
 end
 for i = 1:nx+1
     line([xef(i) xef(i)],[0 1],'Marker','.','LineStyle','--')
 end
 line([ecent(1)-elemdx/2 ecent(1)-elemdx/2],[0 1],'LineWidth',3)
 for i = 1:num_elem
     line([ecent(i)+elemdx/2 ecent(i)+elemdx/2],[0 1],'LineWidth',3)
 end
 hold off
 xlim([0,1])
 ylim([0,1])
 