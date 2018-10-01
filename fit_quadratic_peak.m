function [x, y] = fit_quadratic_peak(img,cx,cy,s)
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
x0 = max(1,cx-s);
y0 = max(1,cy-s);
x1 = min(size(img,2),cx+s);
y1 = min(size(img,1),cy+s);
b = img(y0:y1,x0:x1);

[x,y] = meshgrid(x0:x1,y0:y1);
x = x - cx;
y = y - cy;
% fit the quadratic model 1 + x +y + x^2 + y^2 + x y
x = x(:);
y = y(:);
b = b(:);
A = [ones(size(x)) x y x.^2 y.^2 x.*y];
%X = (A'*A) \ A'*b;
X = A \ b;
% find the maximum (null derivative)
d = -[X(4), X(6);X(6) X(5)] \ [X(2); X(3)];
p = 0.5;
x = cx + p * d(1);
y = cy + p * d(2);

