function b = mydst1(a)
% Copyright 2020 The MathWorks, Inc.

% Apply DST-I in each column
[N1,N2] = size(a);
xx = [zeros(1,N2);
    a;
    zeros(1,N2);
    zeros(N1,N2)]; % DST1
tmp = - fft(xx)*sqrt(2/(N1+1));
b = imag(tmp(2:N1+1,:));