function p = solvePoissonEquation_2dDCT(b,nx,ny,dx,dy)
% Copyright 2020 The MathWorks, Inc.

% modified wavenumber
kx = 0:nx-1;
ky = 0:ny-1;
mwx = 2*(cos(pi*kx/nx)-1)/dx^2;
mwy = 2*(cos(pi*ky/ny)-1)/dy^2;

% 2D DCT of b (Right hand side)
fhat = dct2(b); % Needs Image Processing Toolbox
% fhat = dct(dct(b)')'; % Same as above (Needs Signal Processing Toolbox instead)

[MWX, MWY] = ndgrid(mwx,mwy);
phat = fhat./(MWX+MWY);

% The solution is not unique (uhat(0,0) = inf);
% Here we fix the mean ( with kx=0,ky=0) to be 0
phat(1,1) = 0;

% Inverse 2D DCT
p = idct2(phat); % Needs Image Processing Toolbox
% u = idct(idct(uhat)')'; % Same as above (Needs Signal Processing Toolbox instead)

end
