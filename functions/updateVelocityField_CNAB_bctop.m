function [u,v] = updateVelocityField_CNAB_bctop(u,v,nx,ny,dx,dy,Re,dt,bctop,method)
% Copyright 2020 The MathWorks, Inc.

persistent Nu_old Nv_old

if isempty(Nu_old) || isempty(Nv_old) || any(size(Nu_old) ~= [nx-1,ny])
    Nu_old = zeros(nx-1,ny);
    Nv_old = zeros(nx,ny-1);
end

% Apply boundary conditions:
% represented as the values on ghost cells
u(:,1) = 0-u(:,2); v(:,1) = 0.;             %bottom
u(:,end) = 2*bctop-u(:,end-1);  v(:,end) = 0;  %top
u(1,:) = 0;  v(1,:) = 0-v(2,:);             %left
u(end,:) = 0;  v(end,:) = 0-v(end-1,:);    %right

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
% Implicit treatment for xy direction
Lubc = zeros(size(Luy));
Lubc(:,end) = 2*bctop/dy^2; % effect of the top BC on Lu
Lvbc = zeros(size(Lvy));

b = u(2:end-1,2:end-1) - dt*((3*Nu-Nu_old)/2 - 1/(2*Re)*(Lux+Luy+Lubc));
xu = getIntermediateU_xyCNAB(u, b, dt, Re, nx, ny, dx, dy);
b = v(2:end-1,2:end-1) - dt*((3*Nv-Nv_old)/2 - 1/(2*Re)*(Lvx+Lvy+Lvbc));
xv = getIntermediateV_xyCNAB(v, b, dt, Re, nx, ny, dx, dy);
   
Nu_old = Nu;
Nv_old = Nv;

u(2:end-1,2:end-1) = xu;
v(2:end-1,2:end-1) = xv;

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
        
    otherwise
        error("Specified method: " + method + " is not supported." + ...
            "It should be either direct or dct");
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