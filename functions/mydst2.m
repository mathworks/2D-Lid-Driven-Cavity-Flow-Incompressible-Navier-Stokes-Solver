function b = mydst2(a)
% Copyright 2020 The MathWorks, Inc.

% Apply DST-II in each column
[N1,N2] = size(a);
xx = [a;
    zeros(N1,N2)];

tmp = fft(xx);
xd = tmp(2:N1+1,:);
ww = sqrt(2/N1)*exp(-1i*pi*(1:N1)/(2*N1)).';
b = xd.*ww;
b(end,:) = b(end,:)./sqrt(2);
b = -imag(b);