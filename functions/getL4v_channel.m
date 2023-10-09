% Copyright 2023 The MathWorks, Inc.
function L4v = getL4v_channel(nx,ny,dx,dy,mask)

% Note: Boundary condition effects x derivative only
% x-derivatives
% vb = (v0+v1)/2 => v0 = 2*vb - v1;
% Lv_1 = (v2 - 2*v1 + v0)/dx^2
%      = (v2 - 2*v1 + 2*vb - u1)/dx^2
%      = (v2 - 3*v1)/dx^2 + 2*vb/dx^2
% for peridic case no bc term from x-derivatives
% 
% y-derivatives
% Lv_1 = (v2 - 2*v1 + vb)/dy^2
%      = (v2 - 2*v1)/dy^2 + vb/dy^2

matrixSize = [nx,ny-1];

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
for jj=1:ny-1
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

        tmpx = 2 + sum(~idx(1:2));
        tmpy = 2;
        
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
L4v = sparse(i,j,v,nx*(ny-1),nx*(ny-1));
L4v(:,~any(L4v,1)) = [];
L4v(~any(L4v,2),:) = [];

end