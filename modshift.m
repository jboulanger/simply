function v = modshift(u,k,mode)
% v = modshift(u,k,mode)
%
% Shift modulations
%
% Jerome Boulanger
%

    if mode == 1
        N = size(u,1);
        [x,y] = meshgrid(0:N-1,0:N-1);
        x = x - N / 2 - 1;
        y = y - N / 2 - 1;
        P = exp(1i*2.0*pi/N*(k(2)*x+k(1)*y));
        v = fft2(ifft2(u) .* P);
    elseif mode == 0
        v = circshift(u,round(k));
    else
        N = size(u,1);
        [x,y] = meshgrid(0:N-1,0:N-1);
        v = ifftshift(interp2(fftshift(u),x-k(2),y-k(1),'cubic'));
        v(isnan(v)) = 0;
    end

end
