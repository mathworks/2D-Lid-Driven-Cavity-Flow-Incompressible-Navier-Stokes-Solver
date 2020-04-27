% Copyright 2020 The MathWorks, Inc.
function xu = getIntermediateU_xyRK3_dst(u, b, dt, Re, nx, ny, dx, dy, id)

kx = [1:nx-1]';
ax = pi*kx/(nx);
mwx = 2*(cos(ax)-1)/dx^2;% DST-I

ky = [1:ny]';
ay = pi*ky/(ny);
mwy = 2*(cos(ay)-1)/dy^2; % DST-II

mw = mwx+mwy'; % Modified Wavenumber

% 各 column 毎に変換されるので、
% column 方向（x 方向に）dst1, 転置して
% y 方向に dst2、転置し返して元に戻す処理
rhshat = mydst2(mydst1(b)')';

switch id
    case 1
        uhat = rhshat./(1-(4/15)*dt/(Re)*mw);
    case 2
        uhat = rhshat./(1-(1/15)*dt/(Re)*mw);
    case 3
        uhat = rhshat./(1-(1/6)*dt/(Re)*mw);
end

xu = mydst1(mydst3(uhat')');

end