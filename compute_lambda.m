function [lambda,beta1,beta2,beta3] = compute_lambda(f,uface,scale);

% compute ratio of largest to smallest weno smoothness ratio.
[ny,nf] = size(uface);
[ny2,nx] = size(f);

% Input checking
if ~isequal(ny,ny2); error(['Inputs to weno5_face_input should have ' ... 
                      'same number of rows']); 
end
if ~isequal(nf,nx-5); error(['First input to weno5_face_input should ' ...
                      'have five more columns than the second']); 
end

% compute indices assuming velocity positive everywhere
ixu3 = sub2ind([ny nx],[1:ny]'*ones(1,nf),ones(ny,1)*[1:nf]); 
ixu2 = ixu3 + ny; ixu1 = ixu2 + ny; 
ixd1 = ixu1 + ny; ixd2 = ixd1 + ny; ixd3 = ixd2 + ny;

% find locations of negative velocities and reverse stencil in
% those locations
nind = find(uface<0); ixu3(nind) = ixd3(nind);
tmp = ixu2(nind); ixu2(nind) = ixd2(nind); ixd2(nind) = tmp;
tmp = ixu1(nind); ixu1(nind) = ixd1(nind); ixd1(nind) = tmp;

% compute scaled second derivative everywhere.
secnd = zeros(ny,nx); 
secnd(:,2:end-1) = (13/12)*(f(:,3:end) - 2*f(:,2:end-1) + f(:,1:end-2)).^2;

% Assume positive uface and fill in values of fluxes appropriately.
flxu3 = f(ixu3);   % three points upwind of face.
flxu2 = f(ixu2); % two points upwind of face.
flxu1 = f(ixu1); % one point upwind of face.
flxd1 = f(ixd1); % one point downwind of face.
flxd2 = f(ixd2); % two points downwind of face.

% Estimate the smoothness for each stencil 
beta1 =secnd(ixu2) +0.25*( flxu3 -4*flxu2 +3*flxu1 ).^2;
beta2 =secnd(ixu1) +0.25*( flxu2            -flxd1 ).^2;
beta3 =secnd(ixd1) +0.25*( 3*flxu1 -4*flxd1 +flxd2 ).^2;

lambda = max(beta1, max(beta2,beta3)) ...
        ./(scale^2*1e-8 + min(beta1, min(beta2,beta3)));
