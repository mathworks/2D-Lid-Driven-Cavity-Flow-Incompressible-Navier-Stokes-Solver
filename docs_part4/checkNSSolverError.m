function [L2error, h_fig] = checkNSSolverError(Re,a,N,dt,tEnd,fuFunc,fvFunc,usolFunc,vsolFunc,psolFunc,timeScheme, visFlag, recordGIF)

% Animation Setting
% recordGIF = false; % set to true if you need GIF
% visFlag = false; % set to true to visualize the flow
visRate = 4; % downsample rate of the data for quiver
recordRate = 20;
filename = 'errorCheck.gif'; % Specify the output file name

% Simulation Setting
Lx = 1; Ly = 1; % domain size
Nx = N; Ny = N; % Number of grids

% Grid size (Equispaced)
dx = Lx/Nx;
dy = Ly/Ny;
% Coordinate of each grid (cell center)
xce = ((1:Nx)-0.5)*dx;
yce = ((1:Ny)-0.5)*dy;
% Coordinate of each grid (cell corner)
xco = (0:Nx)*dx;
yco = (0:Ny)*dy;

% cell center mesh grid
[Xce,Yce] = meshgrid(xce,yce); 
%% Initialize the memory for flow field.
u = nan(Nx+1,Ny+2); % velocity in x direction (u)
v = nan(Nx+2,Ny+1); % velocity in y direction (v)
p = nan(Nx+1,Ny+1); % velocity in y direction (v)

% For the current boundary condition
velbc.uTop = nan(size(u,1),1);
velbc.uBottom = nan(size(u,1),1);
velbc.uLeft = nan(1,size(u,2));
velbc.uRight = nan(1,size(u,2));
velbc.vTop = nan(size(v,1),1);
velbc.vBottom = nan(size(v,1),1);
velbc.vLeft = nan(1,size(v,2));
velbc.vRight = nan(1,size(v,2));

% For the boundary conditoin at one time step ahead
velbc1 = velbc;

%% Initialization of the flow field
t = 0;
[uX,uY] = meshgrid(xco(2:end-1),yce);
usol = usolFunc(uX', uY', t, Re, a);

[vX,vY] = meshgrid(xce,yco(2:end-1));
vsol = vsolFunc(vX', vY', t, Re, a);

[pX,pY] = meshgrid(xce,yce);
p = psolFunc(pX', pY', t, Re, a);
u(2:end-1,2:end-1) = usol;
v(2:end-1,2:end-1) = vsol;

% to be consistent with the velocity arrays, u, v
% Top, Buttom: column vector
% Left, Right: row vector
velbc = getVelocityBoudanryCondition(velbc, usolFunc, vsolFunc, xco, yco, xce, yce, t, Re, a);

u(:,1) = 2*velbc1.uBottom - u(:,2); v(:,1) = velbc1.vBottom;             %bottom
u(:,end) = 2*velbc1.uTop - u(:,end-1);  v(:,end) = velbc1.vTop;  %top
u(1,:) = velbc1.uLeft;  v(1,:) = 2*velbc1.vLeft - v(2,:);             %left
u(end,:) = velbc1.uRight;  v(end,:) = 2*velbc1.vRight - v(end-1,:);    %right

% divergence check
inflow = sum(velbc.vBottom(2:end-1))*dx + sum(velbc.uLeft(2:end-1))*dy;
outflow = sum(velbc.vTop(2:end-1))*dx + sum(velbc.uRight(2:end-1))*dy;
assert(abs(inflow - outflow) < eps, "Inflow flux must match the outflow flux.")


if visFlag
   %% Contour plot
    h_fig = figure;
    h_fig.Position = [350  150  700  300];
    ha1 = subplot(1,2,1);
    ha2 = subplot(1,2,2);
    
    % Simulation Results
    % get velocity at the cell center (for visualization)
    uce = (u(1:end-1,2:end-1)+u(2:end,2:end-1))/2; % u at cell center
    vce = (v(2:end-1,1:end-1)+v(2:end-1,2:end))/2; % v at cell center
    [~,h_abs] = contourf(ha1, Xce',Yce',sqrt(uce.^2+vce.^2)); % contour
    
    % Analytical Results
    % get velocity at the cell center (for visualization)
    usol = usolFunc(Xce', Yce', t, Re, a);
    vsol = vsolFunc(Xce', Yce', t, Re, a);
    [~,h_absSol] = contourf(ha2, Xce',Yce',sqrt(usol.^2+vsol.^2)); % contour
    
    hold(ha1,'on');
    hold(ha2,'on');
    
    % Some cosmetics
    h_abs.LevelList = linspace(-1,1,24);
    h_absSol.LevelList = linspace(-1,1,24);
    ha1.CLim = [0,1];
    ha2.CLim = [0,1];
    ha1.FontSize = 14;
    ha2.FontSize = 14;
    ha1.XTick = [];ha1.YTick = [];
    ha2.XTick = [];ha2.YTick = [];
    
    %% Quiver plot
    % Downsample the data（d = downsampled）
    xced = xce(1:visRate:end);
    yced = yce(1:visRate:end);
    [Xced,Yced] = meshgrid(xced, yced);
    
    % Simulation Results
    uced = uce(1:visRate:end,1:visRate:end);
    vced = vce(1:visRate:end,1:visRate:end);
    h_quiver = quiver(ha1, Xced',Yced',uced,vced,1,'Color',[0,0,0]);
    
    % Analytical Results
    usold = usolFunc(Xced', Yced', t, Re, a);
    vsold = vsolFunc(Xced', Yced', t, Re, a);
    h_quiverSol = quiver(ha2, Xced',Yced',usold,vsold,1,'Color',[0,0,0]);
    
    hold(ha1,'off');
    xlim(ha1,[0 Lx]); ylim(ha1,[0 Ly]);
    hold(ha2,'off');
    xlim(ha2,[0 Lx]); ylim(ha2,[0 Ly]);
    
    
    title(ha1,'Simulation Results');
    title(ha2,'Analytical Results');
    axis(ha1,'equal');
    axis(ha2,'equal');
end

%% Start the Simulation
% Just to make it a little fun, the velocity of the top lid reverses at the
% 1000 steps (out of nt steps).

initialFrame = true;
index = 1;
while t < tEnd
    
    velbc = getVelocityBoudanryCondition(velbc1, usolFunc, vsolFunc, xco, yco, xce, yce, t, Re, a);
    velbc1 = getVelocityBoudanryCondition(velbc1, usolFunc, vsolFunc, xco, yco, xce, yce, t+dt, Re, a);
    
    uForce = fuFunc(uX', uY', t, Re, a);
    vForce = fvFunc(vX', vY', t, Re, a);
    %     psol = psolFunc(Xce', Yce', t, Re, a);
    
    % Update the velocity field (uses dct)
    switch timeScheme
        case 'Euler'
            [u,v,p] = updateVelocityField_Euler(u,v,[],Nx,Ny,dx,dy,Re,dt,velbc,'dct',uForce,vForce,velbc1);
        case 'CNAB'
            [u,v,p] = updateVelocityField_CNAB( u,v,Nx,Ny,dx,dy,Re,dt,velbc,'dct',uForce,vForce,velbc1);
        case 'RK3'
            [u,v,p] = updateVelocityField_RK3(u,v,Nx,Ny,dx,dy,Re,dt,velbc,'dct',uForce,vForce,velbc1);
        otherwise
            error('timeScheme is not recognized');
    end
    
    t = t + dt;
    
    % get velocity at the cell center (for visualization)
    uce = (u(1:end-1,2:end-1)+u(2:end,2:end-1))/2; % u at cell center
    vce = (v(2:end-1,1:end-1)+v(2:end-1,2:end))/2; % v at cell center
    usol = usolFunc(Xce', Yce', t, Re, a);
    vsol = vsolFunc(Xce', Yce', t, Re, a);
    
    uSim = [uce(:); vce(:)];
    uSol = [usol(:); vsol(:)];
    L2error = norm(uSim-uSol)/norm(uSol);
    
    CFL = max([u(:)/dx; v(:)/dy])*dt;
%     disp("error: " + num2str(L2error) + " at " + num2str(t) + " CFL: " + num2str(CFL));
    
    % Update the plot at every recordRate steps
    if visFlag && mod(index,recordRate) == 0
        
        % update plot (downsample)
        h_quiver.UData = uce(1:visRate:end,1:visRate:end);
        h_quiver.VData = vce(1:visRate:end,1:visRate:end);
        h_abs.ZData = sqrt(uce.^2+vce.^2);
        
        h_quiverSol.UData = usol(1:visRate:end,1:visRate:end);
        h_quiverSol.VData = vsol(1:visRate:end,1:visRate:end);
        h_absSol.ZData = sqrt(usol.^2+vsol.^2);
        
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
    
    index = index + 1;
    
   
end
 disp("(N,dt) = (" + num2str(N) + ","  + num2str(dt) +") error: " + num2str(L2error) + " at " + num2str(t) + " CFL: " + num2str(CFL));
 
end

function velbc = getVelocityBoudanryCondition(velbc, usolFunc, vsolFunc, xco, yco, xce, yce, t, Re, a)


velbc.uTop = usolFunc(xco, 1, t, Re, a)';
velbc.uBottom = usolFunc(xco, 0, t, Re, a)';
velbc.uLeft(2:end-1) = usolFunc(0, yce, t, Re, a);
velbc.uRight(2:end-1) = usolFunc(1, yce, t, Re, a);
velbc.vTop(2:end-1) = vsolFunc(xce, 1, t, Re, a)';
velbc.vBottom(2:end-1) = vsolFunc(xce, 0, t, Re, a)';
velbc.vLeft = vsolFunc(0, yco, t, Re, a);
velbc.vRight = vsolFunc(1, yco, t, Re, a);

end