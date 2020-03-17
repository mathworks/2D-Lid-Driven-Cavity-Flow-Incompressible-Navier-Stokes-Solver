% Copyright 2020 The MathWorks, Inc.
function Ax = operatorAv_CNAB(x,dt,Re,nx,ny,dx,dy)
% Getting A for v velocity (Crank-Nicolson)
% A = [I - dt/(2*Re)*Lxy];

x = reshape(x,nx,ny-1);
xbig = zeros(nx+2,ny+1);
xbig(2:end-1,2:end-1) = x;

% Addting boudanry value
% ub = (u0+u1)/2 => u0 = 2*ub - u1;
xbig(1,:) = - xbig(2,:);
xbig(end,:) = - xbig(end-1,:);

Ax = x - dt/(2*Re)*(xbig(1:end-2,2:end-1)-2*xbig(2:end-1,2:end-1)+xbig(3:end,2:end-1))/dx^2; % nx-1 * ny
Ax = Ax - dt/(2*Re)*(xbig(2:end-1,1:end-2)-2*xbig(2:end-1,2:end-1)+xbig(2:end-1,3:end))/dy^2; % nx-1 * ny
Ax = Ax(:);
end

