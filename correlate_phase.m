function [xc, yc] = correlate_phase(x,y)
% Phase correlation
%
r = fftn(tape(x)) .* conj(fftn(tape(y)));
r = fftshift(real(ifftn(r./abs(r))));
clf; imshow((log(r)),[])
[yc, xc] = find(r == max(r(:)));
[xc, yc] = fit_quadratic_peak(r,xc,yc,1);
hold on;
plot(xc,yc,'o');
hold off
xc = xc - size(x,2) / 2 - 1;
yc = yc - size(x,1) / 2 - 1;
if nargout == 1
    xc = [xc;yc];
end
