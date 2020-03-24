function b = mydst3(a)
% Copyright 2020 The MathWorks, Inc.

% Apply DST-III in each column
[N1,N2] = size(a);
a(end,:) = a(end,:)./sqrt(2);

xx = [zeros(1,N2);
    a;
    zeros(N1-1,N2)]; % DST1

% weight
ww = exp(-1i*pi*(0:2*N1-1)/(2*N1)).';
xxw = xx.*ww;
tmp = fft(xxw)*sqrt(2/N1);
b = -imag(tmp(1:N1,:));

