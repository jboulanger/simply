function otf = gen_otf(n,p)
  v0 = n/p;
  [x, y] = meshgrid(-n/2 : n/2-1, -n/2 : n/2-1);
  v = min(1,  sqrt(x .* x + y .* y) / v0);
  otf = 2 / pi * (acos(v) - v .* sqrt(1-v.*v));
  otf = fftshift(otf);
  otf = otf ./ max(otf(:));
end
