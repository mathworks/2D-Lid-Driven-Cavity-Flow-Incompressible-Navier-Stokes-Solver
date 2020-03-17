% Copyright 2020 The MathWorks, Inc.
function xv = getIntermediateV_xyCNAB(v, b, dt, Re, nx, ny, dx, dy)

persistent A4v dtRe_old L4v

if isempty(dtRe_old)
    dtRe_old = dt/Re;
end

if isempty(A4v) || any(size(L4v) ~= [nx*(ny-1),nx*(ny-1)]) || dtRe_old ~= dt/Re
    dtRe_old = dt/Re;
    maskV = false(nx+2,ny+1);
    maskV(2:end-1,2:end-1) = true;
    L4v = getL4v(nx,ny,dx,dy,maskV);
    A4v = eye(size(L4v),'like',L4v)-dt/(2*Re)*L4v;
    A4v = decomposition(A4v);  
end

% x0 = v(2:end-1,2:end-1);
% [xv,flag] = cgs(A4v,b(:),[],[],[],[],x0(:));
xv = A4v\b(:);
xv = reshape(xv,[nx,ny-1]);

end

