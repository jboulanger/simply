function [data,otf] = generate_sim_data(obj,m,otf,zoom,noise_std)

L = size(m,3);
data = real(ifft2(fft2(m.*repmat(obj,[1 1 L])) .* repmat(otf,[1 1 L])));
if zoom == 2
    % downsample data;
    data = fftshift(fftshift(fft2(data),1),2);
    data = data(size(data,1)/4+1:3*size(data,1)/4,size(data,2)/4+1:3*size(data,2)/4,:);
    data = real(ifft2(ifftshift(ifftshift(data,2),1)));
    % downsample the otf    
    otf = fftshift(otf);
    otf = otf(size(otf,1)/4+1:3*size(otf,1)/4,size(otf,2)/4+1:3*size(otf,2)/4,:);
    otf = ifftshift(otf);
end
data = data + noise_std * randn(size(data));