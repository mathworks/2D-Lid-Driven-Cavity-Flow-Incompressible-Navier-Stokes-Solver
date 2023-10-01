% Copyright (c) 2023, The MathWorks, Inc.
% This is example code to generate markers from object from image
I = imread('image.jpg'); % An image with object that you want makers to imitate.

BW = imbinarize(im2gray(I),0.7);
BW = imresize(BW, [256,256]);
BW = flipud(BW);

mgrid = linspace(0,1,256);
[mX,mY] = meshgrid(mgrid,mgrid);
mmX = mX(BW);
mmY = mY(BW);

pos = [mmX,mmY];
vel = zeros(size(pos));
markers = table(pos,vel);