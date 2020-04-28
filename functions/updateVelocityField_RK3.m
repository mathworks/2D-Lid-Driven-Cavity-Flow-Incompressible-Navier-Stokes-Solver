% Copyright 2020 The MathWorks, Inc.
function [u,v,p] = updateVelocityField_RK3(u,v,nx,ny,dx,dy,Re,dt,velbc,method,uForce,vForce,velbc1)

[u,v,~] = updateVelocityField_RK3substep(u,v,nx,ny,dx,dy,Re,dt,velbc,method,1,uForce,vForce,velbc1);
[u,v,~] = updateVelocityField_RK3substep(u,v,nx,ny,dx,dy,Re,dt,velbc,method,2,uForce,vForce,velbc1);
[u,v,p] = updateVelocityField_RK3substep(u,v,nx,ny,dx,dy,Re,dt,velbc,method,3,uForce,vForce,velbc1);

% check the divergence
%     b = ((u(2:end,2:end-1)-u(1:end-1,2:end-1))/dx ...
%         + (v(2:end-1,2:end)-v(2:end-1,1:end-1))/dy);
%     norm(b)

end

function [u,v,p] = updateVelocityField_RK3substep(u,v,nx,ny,dx,dy,Re,dt,velbc,method,id,uForce,vForce,velbc1)

persistent Nu_old Nv_old

if isempty(Nu_old) || isempty(Nv_old) || any(size(Nu_old) ~= [nx-1,ny])
    Nu_old = zeros(nx-1,ny);
    Nv_old = zeros(nx,ny-1);
end

% Apply boundary conditions:
% represented as the values on ghost cells
u(:,1) = 2*velbc.uBottom - u(:,2); v(:,1) = velbc.vBottom;             %bottom
u(:,end) = 2*velbc.uTop - u(:,end-1);  v(:,end) = velbc.vTop;  %top
u(1,:) = velbc.uLeft;  v(1,:) = 2*velbc.vLeft - v(2,:);             %left
u(end,:) = velbc.uRight;  v(end,:) = 2*velbc.vRight - v(end-1,:);    %right

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
% forcing
Nu = Nu - uForce;

% 3-2. get derivative for v
Nv = (vvce(:,2:end) - vvce(:,1:end-1))/dy;
Nv = Nv + (uvco(2:end,2:end-1) - uvco(1:end-1,2:end-1))/dx;
% forcing
Nv = Nv - vForce;

%  仮の速度場算出
Lubc = zeros(size(Luy));
Lubc(1,:) = velbc1.uLeft(2:end-1)/dx^2;
Lubc(end,:) = velbc1.uRight(2:end-1)/dx^2;
Lubc(:,1) = 2*velbc1.uBottom(2:end-1)/dy^2;
Lubc(:,end) = 2*velbc1.uTop(2:end-1)/dy^2;

Lvbc = zeros(size(Lvy));
Lvbc(1,:) = 2*velbc1.vLeft(2:end-1)/dx^2;
Lvbc(end,:) = 2*velbc1.vRight(2:end-1)/dx^2;
Lvbc(:,1) = velbc1.vBottom(2:end-1)/dy^2;
Lvbc(:,end) = velbc1.vTop(2:end-1)/dy^2;

alpha = [4/15,1/15,1/6];
% beta = alpha;
gamma = [8/15,5/12,3/4];
zeta = [0,-17/60,-5/12];

% Implicit treatment for xy direction
b = u(2:end-1,2:end-1) - dt*(gamma(id)*Nu + zeta(id)*Nu_old - alpha(id)/(Re)*(Lux+Luy+Lubc));
% xu = getIntermediateU_xyRK3(u, b, dt, Re, nx, ny, dx, dy, id);
xu = getIntermediateU_xyRK3_dst(u, b, dt, Re, nx, ny, dx, dy, id);

b = v(2:end-1,2:end-1) - dt*(gamma(id)*Nv + zeta(id)*Nv_old - alpha(id)/(Re)*(Lvx+Lvy+Lvbc));
% xv = getIntermediateV_xyRK3(v, b, dt, Re, nx, ny, dx, dy, id);
xv = getIntermediateV_xyRK3_dst(v, b, dt, Re, nx, ny, dx, dy, id);

u(2:end-1,2:end-1) = xu;
v(2:end-1,2:end-1) = xv;

Nu_old = Nu;
Nv_old = Nv;


u(:,1) = 2*velbc1.uBottom - u(:,2); v(:,1) = velbc1.vBottom;             %bottom
u(:,end) = 2*velbc1.uTop - u(:,end-1);  v(:,end) = velbc1.vTop;  %top
u(1,:) = velbc1.uLeft;  v(1,:) = 2*velbc1.vLeft - v(2,:);             %left
u(end,:) = velbc1.uRight;  v(end,:) = 2*velbc1.vRight - v(end-1,:);    %right

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
    case 'dct_p'
        % by using the discrete cosine transform（コサイン変換使用）
        % Note: Signal Processing Toolbox required
        dp = solvePoissonEquation_2dDCT_p(b,nx,ny,dx,dy);
        
    otherwise
        error("Specified method: " + method + " is not supported." + ...
            "It should be either direct or dct");
end

% correction to get the final velocity
switch method
    case 'dct_p'
        u(2:end-1,2:end-1) = u(2:end-1,2:end-1) -  (dp(2:end,:)-dp(1:end-1,:))/dx;
        v(2:end-1,2:end-1) = v(2:end-1,2:end-1) -  (dp(:,2:end)-dp(:,1:end-1))/dy;
        u(end,2:end-1) = u(end,2:end-1) + 2*dp(end,:)/dx;    %right
        v(2:end-1,end) = v(2:end-1,end) + 2*dp(:,end)/dy;    %top
        v(2:end-1,1) = v(2:end-1,1) - 2*dp(:,1)/dy;    %bottom
    otherwise
        u(2:end-1,2:end-1) = u(2:end-1,2:end-1) -  (dp(2:end,:)-dp(1:end-1,:))/dx;
        v(2:end-1,2:end-1) = v(2:end-1,2:end-1) -  (dp(:,2:end)-dp(:,1:end-1))/dy;
end

p = dp;

end