addpath('../functions/');
nx = 5;
ny = 5;

maskV = false(nx+2,ny+1);
maskV(1:end,2:end-1) = true;
L4v = getL4v_channel(nx,ny,dx,dy,maskV);
A4v1 = speye(size(L4v))-(4/15)*dt/(Re)*L4v; %beta(1)
A4v2 = speye(size(L4v))-(1/15)*dt/(Re)*L4v; %beta(1)
A4v3 = speye(size(L4v))-(1/6)*dt/(Re)*L4v; %beta(1)
A4v1 = decomposition(A4v1);
A4v2 = decomposition(A4v2);
A4v3 = decomposition(A4v3);
