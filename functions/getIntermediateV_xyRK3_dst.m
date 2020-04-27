% Copyright 2020 The MathWorks, Inc.
function xv = getIntermediateV_xyRK3_dst(v, b, dt, Re, nx, ny, dx, dy,id)


kx = [1:nx]';
ax = pi*kx/(nx);
mwx = 2*(cos(ax)-1)/dx^2;% DST-I

ky = [1:ny-1]';
ay = pi*ky/(ny);
mwy = 2*(cos(ay)-1)/dy^2; % DST-II

mw = mwx+mwy'; % Modified Wavenumber

% 各 column 毎に変換されるので、
% column 方向（x 方向に）dst2, 転置して
% y 方向に dst1、転置し返して元に戻す処理
rhshat = mydst1(mydst2(b)')';

switch id
    case 1
        vhat = rhshat./(1-(4/15)*dt/(Re)*mw);
    case 2
        vhat = rhshat./(1-(1/15)*dt/(Re)*mw);
    case 3
        vhat = rhshat./(1-(1/6)*dt/(Re)*mw);
end

xv = mydst3(mydst1(vhat')');

end

