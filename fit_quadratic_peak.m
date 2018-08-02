function [x, y] = fit_quadratic_peak(img,x0,y0,s)
%
% [x, y, intensity] = fit_quadratic_peak(img,x,y,s)
%
% adjust a quadratic model to a peak in the image at location x,y
% s is the size of the neighborhood on which the fit is performed
%
% [x,y] = meshgrid(1:10,1:10);
% z = exp(-0.5*((x-5.5).^2+(y-5.5).^2)/10);
% [cx,cy] = fit_quadratic_peak(z,5,5,2)
%
%
% Jerome Boulanger 2015

b = img(y0-s:y0+s,x0-s:x0+s);
[x,y] = meshgrid(-s:s,-s:s);
% fit the quadratic model 1 + x +y + x^2 + y^2 + x y
x = x(:);
y = y(:);
b = b(:);
A = [ones(size(x)) x y x.^2 y.^2 x.*y];
X = (A'*A) \ A'*b;
% find the maximum (null derivative)
d = -[X(4), X(6);X(6) X(5)] \ [X(2); X(3)];
p = 0.5;
x = x0 + p * d(1);
y = y0 + p * d(2);

