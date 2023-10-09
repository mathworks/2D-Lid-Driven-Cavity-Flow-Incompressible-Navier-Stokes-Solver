addpath('../functions/');
nx = 5;
ny = 5;

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
