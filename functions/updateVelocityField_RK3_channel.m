% Copyright 2023 The MathWorks, Inc.
function [u,v,p] = updateVelocityField_RK3_channel(u,v,nx,ny,dx,dy,Re,dt,velbc,method,uForce,vForce,velbc1)

[u,v,~] = updateVelocityField_RK3substep(u,v,nx,ny,dx,dy,Re,dt,velbc,method,1,uForce,vForce,velbc1);
[u,v,~] = updateVelocityField_RK3substep(u,v,nx,ny,dx,dy,Re,dt,velbc,method,2,uForce,vForce,velbc1);
[u,v,p] = updateVelocityField_RK3substep(u,v,nx,ny,dx,dy,Re,dt,velbc,method,3,uForce,vForce,velbc1);

% check the divergence
% b = ((u(2:nx+1,2:ny+1)-u(1:nx,2:ny+1))/dx ...
    % + (v(2:nx+1,2:ny+1)-v(2:nx+1,1:ny))/dy);
% norm(b)

end

function [u,v,p] = updateVelocityField_RK3substep(u,v,nx,ny,dx,dy,Re,dt,velbc,method,id,uForce,vForce,velbc1)

persistent Nu_old Nv_old

if isempty(Nu_old) || isempty(Nv_old) || any(size(Nu_old) ~= [nx,ny])
    Nu_old = zeros(nx,ny);
    Nv_old = zeros(nx,ny-1);
end

% Apply boundary conditions:
% represented as the values on ghost cells
% u(:,1) = 2*velbc.uBottom - u(:,2); v(:,1) = velbc.vBottom;             %bottom
% u(:,end) = 2*velbc.uTop - u(:,end-1);  v(:,end) = velbc.vTop;  %top
u(:,1) = 2*velbc.uBottom - u(:,2); v(:,1) = velbc.vBottom;             %bottom
u(:,ny+2) = 2*velbc.uTop - u(:,ny+1);  v(:,ny+1) = velbc.vTop;  %top


% ToDo: need to make it periodic
% u(1,:) = u(end-1,:);  
% u(end,:) = u(2,:);  
% v(1,:) = v(end-1,:);
% v(end,:) = v(2,:); 
u(1,:) = u(nx+1,:);  
u(nx+2,:) = u(2,:);  
v(1,:) = v(nx+1,:);
v(nx+2,:) = v(2,:); 
% u(1,:) = velbc.uLeft;  v(1,:) = 2*velbc.vLeft - v(2,:);             %left
% u(end,:) = velbc.uRight;  v(end,:) = 2*velbc.vRight - v(end-1,:);    %right

%  拡散項(u) Get viscous terms for u
% Lux = (u(1:end-2,2:end-1)-2*u(2:end-1,2:end-1)+u(3:end,2:end-1))/dx^2; % nx * ny
% Luy = (u(2:end-1,1:end-2)-2*u(2:end-1,2:end-1)+u(2:end-1,3:end))/dy^2; % nx * ny
Lux = (u(1:nx,2:ny+1)-2*u(2:nx+1,2:ny+1)+u(3:nx+2,2:ny+1))/dx^2; % nx * ny
Luy = (u(2:nx+1,1:ny)-2*u(2:nx+1,2:ny+1)+u(2:nx+1,3:ny+2))/dy^2; % nx * ny

% 拡散項(v) Get viscous terms for v
% Lvx = (v(1:end-2,2:end-1)-2*v(2:end-1,2:end-1)+v(3:end,2:end-1))/dx^2; % nx * ny-1
% Lvy = (v(2:end-1,1:end-2)-2*v(2:end-1,2:end-1)+v(2:end-1,3:end))/dy^2; % nx * ny-1
Lvx = (v(1:nx,2:ny)-2*v(2:nx+1,2:ny)+v(3:nx+2,2:ny))/dx^2; % nx * ny-1
Lvy = (v(2:nx+1,1:ny-1)-2*v(2:nx+1,2:ny)+v(2:nx+1,3:ny+1))/dy^2; % nx * ny-1

% 対流項の計算
% Get nonlinear terms
% 1. interpolate velocity at cell center/cell cornder
% uce = (u(1:end-1,2:end-1)+u(2:end,2:end-1))/2;
% uco = (u(1:end-1,1:end-1)+u(1:end-1,2:end))/2;
% vco = (v(1:end-1,:)+v(2:end,:))/2;
% vce = (v(2:end-1,1:end-1)+v(2:end-1,2:end))/2;
uce = (u(1:nx+1,2:ny+1)+u(2:nx+2,2:ny+1))/2; % average in x
uco = (u(1:nx+1,1:ny+1)+u(1:nx+1,2:ny+2))/2; % average in y
vco = (v(1:nx+1,:)+v(2:nx+2,:))/2; % average in x
vce = (v(2:nx+2,1:ny)+v(2:nx+2,2:ny+1))/2; % average in y

% 2. multiply
uuce = uce.*uce; % nx+1, ny
uvco = uco.*vco; % nx+1, ny+1
vvce = vce.*vce; % nx+1, ny

% 3-1. get derivative for u
% Nu = (uuce(2:end,:) - uuce(1:end-1,:))/dx;
% Nu = Nu + (uvco(2:end,2:end) - uvco(2:end,1:end-1))/dy;
Nu = (uuce(2:nx+1,:) - uuce(1:nx,:))/dx;
Nu = Nu + (uvco(2:nx+1,2:ny+1) - uvco(2:nx+1,1:ny))/dy;
% forcing
Nu = Nu - uForce;

% 3-2. get derivative for v
% Nv = (vvce(:,2:end) - vvce(:,1:end-1))/dy;
% Nv = Nv + (uvco(2:end,2:end-1) - uvco(1:end-1,2:end-1))/dx;
Nv = (vvce(2:nx+1,2:ny) - vvce(2:nx+1,1:ny-1))/dy;
Nv = Nv + (uvco(2:nx+1,2:ny) - uvco(1:nx,2:ny))/dx;
% forcing
Nv = Nv - vForce;

%  仮の速度場算出
Lubc = zeros(size(Luy)); % nx, ny
% Lubc(1,:) = velbc1.uLeft(2:end-1)/dx^2;
% Lubc(end,:) = velbc1.uRight(2:end-1)/dx^2;
% Lubc(:,1) = 2*velbc1.uBottom(2:end-1)/dy^2;
% Lubc(:,end) = 2*velbc1.uTop(2:end-1)/dy^2;
Lubc(:,1) = 2*velbc1.uBottom(2:nx+1)/dy^2;
Lubc(:,ny) = 2*velbc1.uTop(2:nx+1)/dy^2;

Lvbc = zeros(size(Lvy)); % nx, ny-1
% Lvbc(1,:) = 2*velbc1.vLeft(2:end-1)/dx^2;
% Lvbc(end,:) = 2*velbc1.vRight(2:end-1)/dx^2;
% Lvbc(:,1) = velbc1.vBottom(2:end-1)/dy^2;
% Lvbc(:,end) = velbc1.vTop(2:end-1)/dy^2;
Lvbc(:,1) = velbc1.vBottom(2:nx+1)/dy^2;
Lvbc(:,ny-1) = velbc1.vTop(2:nx+1)/dy^2;

alpha = [4/15,1/15,1/6];
% beta = alpha;
gamma = [8/15,5/12,3/4];
zeta = [0,-17/60,-5/12];

% Implicit treatment for xy direction
b = u(2:end-1,2:end-1) - dt*(gamma(id)*Nu + zeta(id)*Nu_old - alpha(id)/(Re)*(Lux+Luy+Lubc));
xu = getIntermediateU_xyRK3_channel(u, b, dt, Re, nx, ny, dx, dy, id);

b = v(2:end-1,2:end-1) - dt*(gamma(id)*Nv + zeta(id)*Nv_old - alpha(id)/(Re)*(Lvx+Lvy+Lvbc));
xv = getIntermediateV_xyRK3_channel(v, b, dt, Re, nx, ny, dx, dy, id);

% u(2:end-1,2:end-1) = xu;
% v(2:end-1,2:end-1) = xv;
u(2:nx+1,2:ny+1) = xu;
v(2:nx+1,2:ny) = xv;

Nu_old = Nu;
Nv_old = Nv;


u(:,1) = 2*velbc1.uBottom - u(:,2); v(:,1) = velbc1.vBottom;             %bottom
u(:,end) = 2*velbc1.uTop - u(:,end-1);  v(:,end) = velbc1.vTop;  %top

% ToDo: need to make it periodic
% u(1,:) = u(end-1,:);  
% u(end,:) = u(2,:);  
% v(1,:) = v(end-1,:);
% v(end,:) = v(2,:); 
u(1,:) = u(nx+1,:);
u(nx+2,:) = u(2,:);
v(1,:) = v(nx+1,:);
v(nx+2,:) = v(2,:);
% u(1,:) = velbc1.uLeft;  v(1,:) = 2*velbc1.vLeft - v(2,:);             %left
% u(end,:) = velbc1.uRight;  v(end,:) = 2*velbc1.vRight - v(end-1,:);    %right


% 新しい速度場
% 圧力の式（ポワソン方程式）を解いて速度場を質量保存を満たす場に写像。
% velocity correction
% RHS of pressure Poisson eq.
% here's TODO
% b = ((u(2:end,2:end-1)-u(1:end-1,2:end-1))/dx ...
    % + (v(2:end,2:end)-v(2:end,1:end-1))/dy);
b = ((u(2:nx+1,2:ny+1)-u(1:nx,2:ny+1))/dx ...
    + (v(2:nx+1,2:ny+1)-v(2:nx+1,1:ny))/dy);
% norm(b1)

% Solve for p
% For channel direct is the only method available now.
% Could be implement fft based later
switch method
    case 'direct'
        % by directly inverting the matrix（直接法）
        dp = solvePoissonEquation_direct_channel(b,nx,ny,dx,dy);
    otherwise
        error("Specified method: " + method + " is not supported for channel case");
end

% correction to get the final velocity
% u(2:end-1,2:end-1) = u(2:end-1,2:end-1) -  (dp(2:end,:)-dp(1:end-1,:))/dx;
% v(2:end-1,2:end-1) = v(2:end-1,2:end-1) -  (dp(1:end-1,2:end)-dp(1:end-1,1:end-1))/dy;
u(2:nx+1,2:ny+1) = u(2:nx+1,2:ny+1) -  (dp(2:nx+1,:)-dp(1:nx,:))/dx;
v(2:nx+1,2:ny) = v(2:nx+1,2:ny) -  (dp(1:nx,2:ny)-dp(1:nx,1:ny-1))/dy;

p = dp;

% u(1,:) = u(end-1,:);
% u(end,:) = u(2,:);
% v(1,:) = v(end-1,:);
% v(end,:) = v(2,:);
u(1,:) = u(nx+1,:);  
u(nx+2,:) = u(2,:);  
v(1,:) = v(nx+1,:);
v(nx+2,:) = v(2,:); 

% Is this nessesary?
u(:,1) = 2*velbc.uBottom - u(:,2); v(:,1) = velbc.vBottom;             %bottom
u(:,ny+2) = 2*velbc.uTop - u(:,ny+1);  v(:,ny+1) = velbc.vTop;  %top

% Check divergence (only for testing)
% b = ((u(2:nx+1,2:ny+1)-u(1:nx,2:ny+1))/dx ...
%     + (v(2:nx+1,2:ny+1)-v(2:nx+1,1:ny))/dy);
% norm(b)

end