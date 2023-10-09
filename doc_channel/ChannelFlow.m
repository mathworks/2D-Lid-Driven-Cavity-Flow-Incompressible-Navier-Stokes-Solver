%% Animation of the Velocity Field
% Making animation of the flow field is the fun part of CFD (Colorful Fluid
% Dynamics). In addition to the contour plot, let's create an arrow plot that
% represents the velocity using |quiver| function. The code below generates GIF.
% The numerical integration process discussed above is now in the function |updateVelocityField.m.|
%
addpath('../functions/');
clear functions, close all

%% Animation Setting
visRate = 4; % downsample rate of the data for quiver
recordGIF = false; % set to true if you wanna make GIF
recordRate = 2;
filename = 'animation_channel.gif'; % Specify the output file name

% Simulation Setting
timeScheme = 'RK3'; % Euler, CNAB RK3 (RK3 only for now)
Re = 300; % Reynolds number
nt = 5000; % max time steps
Lx = 2; Ly = 1; % domain size
Nx = 64; Ny = 32; % Number of grids
dt = 0.01; % time step;
% Grid size (Equispaced)
dx = Lx/Nx;
dy = Ly/Ny;
% Coordinate of each grid (cell center)
xce = ((1:Nx)-0.5)*dx;
yce = ((1:Ny)-0.5)*dy;
% Coordinate of each grid (cell corner)
xco = (0:Nx)*dx;
yco = (0:Ny)*dy;


%% Initialize the flow field.
u = zeros(Nx+2,Ny+2); % velocity in x direction (u)
v = zeros(Nx+2,Ny+1); % velocity in y direction (v)
uce = (u(1:end-2,2:end-1)+u(2:end-1,2:end-1))/2; % u at cell center
vce = (v(2:end-1,1:end-1)+v(2:end-1,2:end))/2; % v at cell center
% Setting up Visualization
% Downsample the velocity field at the interval of |visRate| so that the plot
% is not going to be filled with arrows. Also GIF is updated (appended) at every
% |recordRate| time steps.

velbc.uTop = nan(size(u,1),1);
velbc.uBottom = nan(size(u,1),1);
velbc.uLeft = nan(1,size(u,2));
velbc.uRight = nan(1,size(u,2));
velbc.vTop = nan(size(v,1),1);
velbc.vBottom = nan(size(v,1),1);
velbc.vLeft = nan(1,size(v,2));
velbc.vRight = nan(1,size(v,2));

velbc1 = velbc;
%%

% to be consistent with the velocity arrays, u, v
% Top, Buttom: column vector
% Left, Right: row vector
velbc.uTop(:) = 0;
velbc.uBottom(:) = 0;
velbc.vTop(2:end-1) = 0;
velbc.vBottom(2:end-1) = 0;

% No need: Periodic in x
% velbc.uLeft(:) = 1;
% velbc.uRight(:) = 0;
% velbc.vLeft(:) = 0;
% velbc.vRight(:) = 0;

%% Contour plot

figure;
% Initial velocity field.
u(2:end-1,2:end-1) = 0; v(2:end-1,2:end-1)  = 0;

% get velocity at the cell center (for visualization)
uce = (u(1:end-2,2:end-1)+u(2:end-1,2:end-1))/2; % u at cell center
vce = (v(2:end-1,1:end-1)+v(2:end-1,2:end))/2; % v at cell center
[Xce,Yce] = meshgrid(xce,yce); % cell center grid
[~,h_abs] = contourf(Xce',Yce',uce); % contour

hold on
%% Quiver plot
% Downsample the data（d = downsampled）
% xced = xce(1:visRate:end);
% yced = yce(1:visRate:end);
% [Xced,Yced] = meshgrid(xced, yced);
% uced = uce(1:visRate:end,1:visRate:end);
% vced = vce(1:visRate:end,1:visRate:end);
% h_quiver = quiver(Xced',Yced',uced,vced,1,'Color',[1,1,1]);
% h_surf = surf(u);

hold off
xlim([0 Lx]); ylim([0 Ly]);

% Delete the ticks/tick labels.
haxes = gca;
haxes.XTick = [];
haxes.YTick = [];


%% Start the Simulation
initialFrame = true;

t = 0;
CFL = max([u(:)/dx; v(:)/dy])*dt;
for ii = 1:nt

    if CFL > 1.2
        % // Control for cfl.
        omega_cfl = 0.7; % Under-relaxation
        want_cfl  = 1.0; % Desired Courant number
        dt = dt*(1.0 + omega_cfl*( want_cfl/CFL - 1.0 ));
    end


    % Flow is driven by pressure gradient in x direction
    uForce = ones(Nx,Ny);
    % uForce = zeros(Nx,Ny);
    vForce = zeros(Nx,Ny-1);

    velbc1 = velbc;
    % Update the velocity field (uses dct)
    % Periodic scheme is only implemented in RK3 
    % Euler, CNAB is not yet ready.
    switch timeScheme
        case 'RK3'
            [u,v,p] = updateVelocityField_RK3_channel(u,v,Nx,Ny,dx,dy,Re,dt,velbc,'direct',uForce,vForce,velbc1);
        otherwise
            error('Specified timeScheme is not recognized/supported');
    end

    velbc1.vTop = v(:,1);
    velbc1.vBottom = v(:,end);

    t = t + dt;
    
    % Update the plot at every recordRate steps
    if mod(ii,recordRate) == 0
        % get velocity at the cell center (for visualization)
        uce = (u(1:end-2,2:end-1)+u(2:end-1,2:end-1))/2; % u at cell center
        vce = (v(2:end-1,1:end-1)+v(2:end-1,2:end))/2; % v at cell center
        % update plot (downsample)
        % h_quiver.UData = uce(1:visRate:end,1:visRate:end);
        % h_quiver.VData = vce(1:visRate:end,1:visRate:end);
        % h_abs.ZData = sqrt(uce.^2+vce.^2);
        h_abs.ZData = uce;

        CFL = max([u(:)/dx; v(:)/dy])*dt;
        disp(" at " + num2str(t) + " CFL: " + num2str(CFL));
        
        drawnow
        
        if recordGIF
            frame = getframe(gcf); %#ok<UNRCH> % Figure 画面をムービーフレーム（構造体）としてキャプチャ
            tmp = frame2im(frame); % 画像に変更
            [A,map] = rgb2ind(tmp,256); % RGB -> インデックス画像に
            if initialFrame
                imwrite(A,map,filename,'gif','LoopCount',Inf,'DelayTime',0.1);
                initialFrame = false;
            else
                imwrite(A,map,filename,'gif','WriteMode','append','DelayTime',0.1);% 画像をアペンド
            end
        end 
    end
end

