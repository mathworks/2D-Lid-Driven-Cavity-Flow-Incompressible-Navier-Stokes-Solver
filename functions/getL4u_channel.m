% Copyright 2023 The MathWorks, Inc.
function L4u = getL4u_channel(nx,ny,dx,dy,mask)

% Note: Boundary condition effects y derivative only
% y-derivatives
% ub = (u0+u1)/2 => u0 = 2*ub - u1;
% Lu_1 = (u2 - 2*u1 + u0)/dy^2
%      = (u2 - 2*u1 + 2*ub - u1)/dy^2
%      = (u2 - 3*u1)/dx^2 + 2*ub/dy^2
% 
% x-derivatives
% Lu_1 = (u2 - 2*u1 + ub)/dx^2
%      = (u2 - 2*u1 + u(end))/dx^2 (+ ub/dx^2)
% where u(end) = u(2), u(1) = u(end-1)

matrixSize = [nx,ny];

% 3
%152
% 4
coefx = 1/dx^2;
coefy = 1/dy^2;

i = zeros(nx*ny*5,1);
j = zeros(nx*ny*5,1);
v = zeros(nx*ny*5,1);
index1 = 1;
index = 1;
for jj=1:ny
    for ii=1:nx
        idx = [mask(ii+2,jj+1),mask(ii,jj+1),mask(ii+1,jj+2),mask(ii+1,jj),mask(ii+1,jj+1)];
        stencils = [ii+1,jj
            ii-1,jj
            ii,jj+1
            ii,jj-1
            ii,jj];
        stencils = stencils(idx,:);
        
        % periodic
        x_stencils = stencils(:,1);
        x_stencils = x_stencils - nx*(x_stencils>nx);
        x_stencils(x_stencils==0) = nx;
        stencils(:,1) = x_stencils;

        tmpx = 2;
        tmpy = 2 + sum(~idx(3:4));
        
        coeffs = [coefx
            coefx
            coefy
            coefy
            -tmpx*coefx-tmpy*coefy];
        coeffs = coeffs(idx);
        linearIdx = sub2ind(matrixSize, stencils(:,1), stencils(:,2));
        
        i(index1:index1+length(linearIdx)-1) = index;
        j(index1:index1+length(linearIdx)-1) = linearIdx;
        v(index1:index1+length(linearIdx)-1) = coeffs;
        index1 = index1 + length(linearIdx);
        index = index + 1;
    end
end

i(i==0) = [];
j(j==0) = [];
v(v==0) = [];
L4u = sparse(i,j,v,(nx)*ny,(nx)*ny);
L4u(:,~any(L4u,1)) = [];
L4u(~any(L4u,2),:) = [];

end
