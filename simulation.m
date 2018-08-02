% simulation of structured illumination microscopy
clear all
warning off

zoom = 2;

obj = double(imread('Obj.tif'));
Na = 3;
Np = 3;
cutoff = zoom*2.35;
otf = gen_otf(size(obj,1), cutoff);
period = zoom* 3.16 * ones(Na,Np);
theta = kron((0:Na-1)/Na*180, ones(Np,1));
shift = kron(ones(1,Na), (0:Np-1)'/Np);
L = numel(period);
%
noise_std = 5;


m = double(0.5+0.5*genpat(size(obj,2),size(obj,1),period,theta,shift,1));
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
subplot(121), imshow(mean(data,3),[])

wiener_factor = .001;
padding_level = 0;
prefilter = false;
notchfilter = 1;
subpixel_localization = true;
verbose_mode = 2;

tic
[img,dec,avg] = simrec(data,otf,zoom,wiener_factor,padding_level,prefilter,notchfilter,subpixel_localization,verbose_mode);
toc

figure(1);
otf = fftupscale(otf);
subplot(231), imshow(sqrt(obj),[]);title('WF');
subplot(234), fftshow(obj,otf);title('|FFT|');
subplot(232), imshow(sqrt(dec),[]);title('Deconvolution');
subplot(235), fftshow(dec,otf);title('|FFT|');
subplot(233), imshow(sqrt(img),[]);title('SIM');
subplot(236), fftshow(img,otf);title('|FFT|');

fprintf(1,'dec rms:%f\n', sqrt(mean((dec(:)-obj(:)).^2)));
fprintf(1,'sim rms:%f\n', sqrt(mean((img(:)-obj(:)).^2)));
figure(2),clf,imshow(max(0,img),[])