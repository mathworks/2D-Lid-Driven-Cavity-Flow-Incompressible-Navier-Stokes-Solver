addpath('../functions/');

nx = 5; ny = 2;

tmp = -2*ones(nx,1);
% tmp([1,end]) = -1; % Neumann 境界条件を考慮し両端は -1 (not -2)
% Periodic 境界条件はすべて -2
Ad = diag(tmp);
Au = diag(ones(nx-1,1),1);
Al = diag(ones(nx-1,1),-1);
Ax = Ad+Au+Al;
% for periodic case
Ax(end,1) = 1;
Ax(1,end) = 1;
% Ax(:,:) = 0;

% 上と同じ行列をもう少し完結に書くとこちら。
% まずブロック対角行列成分を作成
dd = eye(nx);
tmp = repmat({sparse(Ax/dx^2 - 2*dd/dy^2)},ny,1);
tmp{1} = Ax/dx^2 - dd/dy^2;
tmp{ny} = Ax/dx^2 - dd/dy^2; % y-方向の端は注意（Neumann BC)
Abig = blkdiag(tmp{:});

% 1ブロック分小さい対角行列で y微分成分を作成
d4y = eye(nx*(ny-1),'like',sparse(1));
Abig(1:end-nx,nx+1:end) = Abig(1:end-nx,nx+1:end) + d4y/dy^2; % 上側
Abig(nx+1:end,1:end-nx) = Abig(nx+1:end,1:end-nx) + d4y/dy^2; % 下側
%Velocity and pressure correction using phi