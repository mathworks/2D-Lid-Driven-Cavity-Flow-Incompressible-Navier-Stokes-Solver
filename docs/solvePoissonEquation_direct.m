function dp = solvePoissonEquation_direct(b,nx,ny,dx,dy)

persistent Abig

if isempty(Abig) || any(size(Abig) ~= [nx*ny,nx*ny])
    
    % 行列 A の構築（メモリ節約のためスパース行列で定義）
    % まず x 方向の微分に関する部分（三重対角行列）から定義します。
    tmp = -2*ones(nx,1);
    tmp([1,end]) = -1; % Neumann 境界条件を考慮し両端は -1 (not -2)
    Ad = diag(tmp);
    Au = diag(ones(nx-1,1),1);
    Al = diag(ones(nx-1,1),-1);
    Ax = Ad+Au+Al;
    
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
    
    % Abig は特異行列でありこのままでは解が一意に定まらないので、
    % 1点を u = 0 と固定して解とします。
    Abig(1,:) = 0;
    Abig(1,1) = 1;
    
end

% 右辺
f = b(:); % ベクトル化

% Abig は特異行列でありこのままでは解が一意に定まらないので、
% 1点を u = 0 と固定して解とします。
f(1) = 0;

% 求解
dp = Abig\f;

% 2次元配列に戻しておきます。
dp = reshape(dp,[nx,ny]);

end