% Copyright 2023 The MathWorks, Inc.
function xu = getIntermediateU_xyRK3_channel(u, b, dt, Re, nx, ny, dx, dy, id)

persistent A4u1 A4u2 A4u3 dtRe_old L4u

if isempty(dtRe_old)
    dtRe_old = dt/Re;
end

if isempty(A4u1) || any(size(L4u) ~= [nx*ny,nx*ny]) || dtRe_old ~= dt/Re
    dtRe_old = dt/Re;
    maskU = false(nx+2,ny+2);
    maskU(1:end,2:end-1) = true;
    L4u = getL4u_channel(nx,ny,dx,dy,maskU);
    A4u1 = speye(size(L4u))-(4/15)*dt/(Re)*L4u; %beta(1)
    A4u2 = speye(size(L4u))-(1/15)*dt/(Re)*L4u; %beta(2)
    A4u3 = speye(size(L4u))-(1/6)*dt/(Re)*L4u; %beta(3)
    A4u1 = decomposition(A4u1);
    A4u2 = decomposition(A4u2);
    A4u3 = decomposition(A4u3);
    
end

% x0 = u(2:end-1,2:end-1);
% [xu,flag] = cgs(A4u,b(:),[],[],[],[],x0(:));
switch id
    case 1
        xu = A4u1\b(:);
    case 2
        xu = A4u2\b(:);
    case 3
        xu = A4u3\b(:);
end

xu = reshape(xu,[nx,ny]);

end