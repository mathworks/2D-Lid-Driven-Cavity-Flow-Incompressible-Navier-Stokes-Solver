% Copyright 2020 The MathWorks, Inc.
function xu = getIntermediateU_xyCNAB(u, b, dt, Re, nx, ny, dx, dy)

persistent A4u dtRe_old L4u

if isempty(dtRe_old)
    dtRe_old = dt/Re;
end

if isempty(A4u) || any(size(L4u) ~= [(nx-1)*ny,(nx-1)*ny]) || dtRe_old ~= dt/Re
    dtRe_old = dt/Re;
    maskU = false(nx+1,ny+2);
    maskU(2:end-1,2:end-1) = true;
    L4u = getL4u(nx,ny,dx,dy,maskU);
    A4u = eye(size(L4u),'like',L4u)-dt/(2*Re)*L4u;
    A4u = decomposition(A4u); 
end

% x0 = u(2:end-1,2:end-1);
% [xu,flag] = cgs(A4u,b(:),[],[],[],[],x0(:));
xu = A4u\b(:);
xu = reshape(xu,[nx-1,ny]);

end