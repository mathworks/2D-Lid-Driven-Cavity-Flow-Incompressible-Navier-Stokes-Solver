% Copyright 2020 The MathWorks, Inc.
function xv = getIntermediateV_xyRK3(v, b, dt, Re, nx, ny, dx, dy,id)

persistent A4v1 A4v2 A4v3 dtRe_old L4v

if isempty(dtRe_old)
    dtRe_old = dt/Re;
end

if isempty(A4v1) || any(size(L4v) ~= [nx*(ny-1),nx*(ny-1)]) || dtRe_old ~= dt/Re
    dtRe_old = dt/Re;
    maskV = false(nx+2,ny+1);
    maskV(2:end-1,2:end-1) = true;
    L4v = getL4v(nx,ny,dx,dy,maskV);
    A4v1 = eye(size(L4v),'like',L4v)-(4/15)*dt/(Re)*L4v; %beta(1)
    A4v2 = eye(size(L4v),'like',L4v)-(1/15)*dt/(Re)*L4v; %beta(1)
    A4v3 = eye(size(L4v),'like',L4v)-(1/6)*dt/(Re)*L4v; %beta(1)
    A4v1 = decomposition(A4v1);
    A4v2 = decomposition(A4v2);
    A4v3 = decomposition(A4v3);
    
end

switch id
    case 1
        xv = A4v1\b(:);
    case 2
        xv = A4v2\b(:);
    case 3
        xv = A4v3\b(:);
end
xv = reshape(xv,[nx,ny-1]);

end

