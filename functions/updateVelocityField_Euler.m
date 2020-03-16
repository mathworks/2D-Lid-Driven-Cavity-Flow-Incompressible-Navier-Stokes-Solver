function [u,v] = updateVelocityField_Euler(u,v,nx,ny,dx,dy,Re,dt,bctop,method)
% Copyright (c) 2020, The MathWorks, Inc.

% Apply boundary conditions:
% represented as the values on ghost cells
u(:,1) = -u(:,2); v(:,1) = 0.;             %bottom
u(:,end) = 2*bctop-u(:,end-1);  v(:,end) = 0.;  %top
u(1,:) = 0.;  v(1,:) = -v(2,:);             %left
u(end,:) = 0.;  v(end,:) = -v(end-1,:);    %right

%  拡散項(u) Get viscous terms for u
Lux = (u(1:end-2,2:end-1)-2*u(2:end-1,2:end-1)+u(3:end,2:end-1))/dx^2; % nx-1 * ny
Luy = (u(2:end-1,1:end-2)-2*u(2:end-1,2:end-1)+u(2:end-1,3:end))/dy^2; % nx-1 * ny

% 拡散項(v) Get viscous terms for v
Lvx = (v(1:end-2,2:end-1)-2*v(2:end-1,2:end-1)+v(3:end,2:end-1))/dx^2; % nx * ny-1
Lvy = (v(2:end-1,1:end-2)-2*v(2:end-1,2:end-1)+v(2:end-1,3:end))/dy^2; % nx * ny-1

% 対流項の計算
% Get nonlinear terms
% 1. interpolate velocity at cell center/cell cornder
uce = (u(1:end-1,2:end-1)+u(2:end,2:end-1))/2;
uco = (u(:,1:end-1)+u(:,2:end))/2;
vco = (v(1:end-1,:)+v(2:end,:))/2;
vce = (v(2:end-1,1:end-1)+v(2:end-1,2:end))/2;

% 2. multiply
uuce = uce.*uce;
uvco = uco.*vco;
vvce = vce.*vce;

% 3-1. get derivative for u
Nu = (uuce(2:end,:) - uuce(1:end-1,:))/dx;
Nu = Nu + (uvco(2:end-1,2:end) - uvco(2:end-1,1:end-1))/dy;

% 3-2. get derivative for v
Nv = (vvce(:,2:end) - vvce(:,1:end-1))/dy;
Nv = Nv + (uvco(2:end,2:end-1) - uvco(1:end-1,2:end-1))/dx;

%  仮の速度場算出
% 一次精度のオイラー積分
% Get intermidiate velocity
u(2:end-1,2:end-1) = u(2:end-1,2:end-1) + dt*(-Nu + (Lux+Luy)/Re);
v(2:end-1,2:end-1) = v(2:end-1,2:end-1) + dt*(-Nv + (Lvx+Lvy)/Re);

% 新しい速度場
% 圧力の式（ポワソン方程式）を解いて速度場を質量保存を満たす場に写像。
% velocity correction
% RHS of pressure Poisson eq.
b = ((u(2:end,2:end-1)-u(1:end-1,2:end-1))/dx ...
    + (v(2:end-1,2:end)-v(2:end-1,1:end-1))/dy);

% Solve for p
switch method
    case 'direct'
        % by directly inverting the matrix（直接法）
        dp = solvePoissonEquation_direct(b,nx,ny,dx,dy);
    case 'dct'
        % by using the discrete cosine transform（コサイン変換使用）
        % Note: Signal Processing Toolbox required
        dp = solvePoissonEquation_2dDCT(b,nx,ny,dx,dy);
    case 'iterative'
        [dp,~] = minres(@(x) operatorDG(x,nx,ny,dx,dy),b(:),1e-4,300);
        dp = reshape(dp,nx,ny);
        
    otherwise
        error("Specified method: " + method + " is not supported." + ...
            "It should be either direct, dct, or iterative");
end
% correction to get the final velocity
p = dp;
u(2:end-1,2:end-1) = u(2:end-1,2:end-1) -  (p(2:end,:)-p(1:end-1,:))/dx;
v(2:end-1,2:end-1) = v(2:end-1,2:end-1) -  (p(:,2:end)-p(:,1:end-1))/dy;

% check the divergence
%     b = ((u(2:end,2:end-1)-u(1:end-1,2:end-1))/dx ...
%         + (v(2:end-1,2:end)-v(2:end-1,1:end-1))/dy);
%     norm(b)

end


function DGx = operatorDG(x,nx,ny,dx,dy)

x = reshape(x,nx,ny);
xbig = zeros(nx+2,ny+2);
xbig(2:end-1,2:end-1) = x;
xbig(1,:) = xbig(2,:);
xbig(end,:) = xbig(end-1,:);
xbig(:,1) = xbig(:,2);
xbig(:,end) = xbig(:,end-1);

DGx = (xbig(1:end-2,2:end-1)-2*xbig(2:end-1,2:end-1)+xbig(3:end,2:end-1))/dx^2; % nx * ny
DGx = DGx + (xbig(2:end-1,1:end-2)-2*xbig(2:end-1,2:end-1)+xbig(2:end-1,3:end))/dy^2; % nx * ny
DGx = DGx(:);

end