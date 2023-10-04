% Copyright (c) 2023 The MathWorks, Inc.
clear
close all
addpath('../functions/');

showContour = true;
recordGIF = false;
recordRate = 30;
filename = 'sampleMarker.gif'; % Specify the output file name

% Customized colormap for vorticity contour
load colormap.mat
% For vorticity contour, good value varies.
maxLevel = 6; % 6 is good for Re=1e4;

Re = 10000; % Reynolds number
nt = 50; % max time steps
Lx = 2; Ly = 1; % domain size
Nx = 160; Ny = 80; % Number of grids
dt = 0.02; % time step;
% ドメイン設定項目はここまで。

%グリッドサイズやグリッドの中央（xce/yce）とコーナー（xco/yco）に当たる部分の座標位置を計算します。
% Grid size (Equispaced)
dx = Lx/Nx;
dy = Ly/Ny;
% Coordinate of each grid (cell center)
xce = ((1:Nx)-0.5)*dx;
yce = ((1:Ny)-0.5)*dy;
% Coordinate of each grid (cell corner)
xco = (0:Nx)*dx;
yco = (0:Ny)*dy;

% Initial condition
u = zeros(Nx+1,Ny+2); % velocity in x direction (u)
v = zeros(Nx+2,Ny+1); % velocity in y direction (v)
p = zeros(Nx,Ny); % pressure (lagurange multiplier)

uce = (u(1:end-1,2:end-1)+u(2:end,2:end-1))/2; % u at cell center
vce = (v(2:end-1,1:end-1)+v(2:end-1,2:end))/2; % v at cell center

% Figure setup
f = figure(1);
[Xce,Yce] = ndgrid(xce,yce); % cell centerの座標グリッド
[Xco,Yco] = ndgrid(xco,yco); % cell cornerの座標グリッド

% One sample marker for setup
pos = [0.1,0.1]; % position
% pos = rand(500,2); % could use this for fun
vel = zeros(size(pos)); % marker velocity (0)
markers = table(pos,vel); % put them in table;

% Passive markers by scatter
hmarker = scatter(markers.pos(:,1),markers.pos(:,2),1,'MarkerFaceColor','red');
hold on
xlim([0 Lx]); ylim([0 Ly]);
grid off

% cosmetic changes on figure
haxes = gca;
haxes.XTick = [];
haxes.YTick = [];
% Figure size
f.Position = [400,300,640,320];
haxes.Position = [0.05,0.05,0.9,0.9];

% Boundary conditions
velbc.uTop = zeros(size(u,1),1);
velbc.uBottom = zeros(size(u,1),1);
velbc.uLeft = zeros(1,size(u,2));
velbc.uRight = zeros(1,size(u,2));
velbc.vTop = zeros(size(v,1),1);
velbc.vBottom = zeros(size(v,1),1);
velbc.vLeft = zeros(1,size(v,2));
velbc.vRight = zeros(1,size(v,2));
velbc1 = velbc;


%% Main simulation begins here
% Adding Forcing term in x direction in a small region
% This is the main driver of the flow.
uForce = zeros(Nx-1,Ny);
vForce = zeros(Nx,Ny-1);

Ny2 = round(Ny/2);
uForce(2:7,Ny2-4:Ny2+5) = 1;

for ii = 1:5000
    % adding Random perturbation at x = 5, y = Ny2;
    vForce(5,Ny2) = (rand-0.5)*0.01;

    % Update velocity field with 3 step Runge-Kutta
    [u,v,p] = updateVelocityField_RK3( u,v,Nx,Ny,dx,dy,Re,dt,velbc,'dct',uForce,vForce,velbc1);

    % Velocity at the cell center
    uce = (u(1:end-1,2:end-1)+u(2:end,2:end-1))/2; % u at cell center
    vce = (v(2:end-1,1:end-1)+v(2:end-1,2:end))/2; % v at cell center

    % Velocity at the marker position (interpolation)
    markerU = interp2(Xce',Yce',uce',markers.pos(:,1),markers.pos(:,2));
    markerV = interp2(Xce',Yce',vce',markers.pos(:,1),markers.pos(:,2));

    % update marker pos (Euler method): rough!
    markers.pos = markers.pos + [markerU, markerV]*dt;

    % Adding more markers!!
    N = 20;
    ymark = linspace(0.47,0.53,N)';
    pos = [0.05*ones(N,1),ymark]; % adding N markers
    vel = zeros(size(pos)); % of velocity 0
    markers = [markers; table(pos,vel)];


    % for recroding animation
    if mod(ii-1,recordRate) == 0

        if showContour
            % Vorticity = dv/dx - du/dy;
            vorticity = (v(2:end,:)-v(1:end-1,:))/dx - (u(:,2:end)-u(:,1:end-1))/dy;

            if ii == 1
                % # levels of contour = 150
                [~,hcontour] = contourf(Xco,Yco,vorticity,150,"LineStyle",'none');
                hcontour.FaceAlpha = 0.75;
                clim([-maxLevel,maxLevel])
                colormap(CustomColormap)

                % make contour plot back of the particles
                h_axes = gca;
                h_axes.Children = h_axes.Children([2,1]);
            else
                hcontour.ZData = vorticity;
                maxV = max(vorticity,[],"all");
                minV = min(vorticity,[],"all");
                % reset levels of contour of 150
                hcontour.LevelList = linspace(minV,maxV,151);
            end

        end

        hmarker.XData = markers.pos(:,1);
        hmarker.YData = markers.pos(:,2);
        drawnow

        if recordGIF
            frame = getframe(gcf); %#ok<UNRCH> % Figure 画面をムービーフレーム（構造体）としてキャプチャ
            tmp = frame2im(frame); % 画像に変更
            % tmp = imresize(tmp,0.5); % To reduce the file size
            [A,map] = rgb2ind(tmp,256); % RGB -> index image

            if ii==1
                imwrite(A,map,filename,'gif','LoopCount',Inf,'DelayTime',0.05);
            else
                imwrite(A,map,filename,'gif','WriteMode','append','DelayTime',0.05);% 画像をアペンド
            end
        end
    end
end
