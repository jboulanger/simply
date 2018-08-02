function [su,dec,avg] = simrec(data,otf,zoom,wiener_factor,padding_level,prefilter,notchfilter,subpixel_localization,verbose_mode)
% [sim,dec,avg] = SIM3X3REC2(data,otf)
% 2D structured illumunation microscopy reconstruction
%
% Jerome Boulanger 2017

if nargin < 3
    zoom = 2;
end
if nargin < 4
    wiener_factor = 0.5;
end
if nargin < 5
    padding_level = 1;
end
if nargin < 6
    prefilter = false;
end
if nargin < 7
    notchfilter = false;
end
if nargin < 8
    subpixel_localization = true;
end
if nargin < 9
    verbose_mode = 2;
end

corrint = 0;
method = 'gust';
L = size(data,3);

% Padding size on each side
padsize = ((2^(9 + padding_level)) - 2^9)/2; 

data = double(data);

tmp = mean(data,3);
data_std = max(tmp(:)) - min(tmp(:));
noise_std = noise_level(tmp);


% Find the order of the images from the peaks in the |fft|
[data,Na,Np] = reorderimages(data, otf, verbose_mode);
%idx = reshape((1:25),5,5)';
%data = data(:,:,idx(:)');
%Na = 5;
%Np = 5;

fprintf('SIM %d angles %d phases\n', Na, Np);

if corrint == 1
    disp('Correct variations of intensities');
    for k = 1:Na
        X = mean(data(:,:,Np*(k-1)+1:Np*k),3);
        for l = 1 : Np
            Y = data(:,:,Np*(k-1)+l);
            p = polyfit(X(:), Y(:), 1);
            data(:,:,Np*(k-1)+l) = (Y - p(2)) / p(1);
        end
    end
end

% Apply a window to prevent edge artifacts
for l = 1:size(data,3)
   data(:,:,l) = tape(data(:,:,l),0.1);
end


% Compute an apotome image
%apo = apotome(data);

if verbose_mode  > 0
    fprintf(1, 'FFT zero padding level %d / %dx%d\n', ...
        padding_level, 2*padsize+2^9, 2*padsize+2^9);
    fprintf(1, 'Noise level %.2f\n', noise_std);
end

if padsize > 0
    data = padarray(data,[padsize padsize 0]);
    otf = abs(ifft2(ifftshift(padarray(fftshift(fft2(otf)),[padsize padsize]))));   
end

otf = otf ./ max(otf(:));
[k,l] = find(otf == 1);
otf = circshift(otf,[k,l]-1);

N = size(data,1);

cutoff = find(otf(1,1:N/2)< eps, 1);
if isempty(cutoff)
    cutoff = N / 2;
end

su = zeros(zoom * N, zoom * N, L);
sw = zeros(zoom * N, zoom * N, L);

if (zoom == 2)
    H = fftupscale(otf);    
else
    H = otf;
end

% mask to detect peaks in the fft
D = H > 1e-3 & H < 0.6;

% wiener parameter = 1/SNR
v = wiener_factor;% * (noise_std / data_std);
if prefilter
    W = conj(H) ./ (v + abs(H).^2);
end
fprintf('wiener parameter %g (SNR = %f / %f = %fdB)\n', v, data_std, noise_std, 20 * log10(data_std / noise_std));
notch =  weights(zoom*N,[0,0],notchfilter,8);

for k = 1:Na
    
    idx = Np*(k-1)+1:Np*k;
    S =  fftn(data(:,:,idx));
    %S = fft2(phasesplit(data(:,:,idx)));
    if zoom == 2
        S = fftupscale(S);
    end
    
    for l = 1:Np
        % pre-wiener
        if prefilter
            S(:,:,l) = W .* S(:,:,l);
        end
        if l == 1
            su(:,:,Np*(k-1)+l) = S(:,:,l);
            sw(:,:,Np*(k-1)+l) = H;
        else         
            % estimate the shift by fitting a peak in fourier space
            %f = peak(abs(S(:,:,l)) .* D, subpixel_localization);
            f = peak2(S(:,:,l),S(:,:,1),D,subpixel_localization);            
            % shift and upscale the components           
            su(:,:,Np*(k-1)+l) = modshift(S(:,:,l).*notch, -f, 1);
            sw(:,:,Np*(k-1)+l) = modshift(H.*notch, -f, 0);
            % correct flip of image
            su(:,:,Np*(k-1)+l) = register_int(sw(:,:,Np*(k-1)+1),sw(:,:,Np*(k-1)+l),su(:,:,Np*(k-1)+1),su(:,:,Np*(k-1)+l));            
            % print information on the components
            period = N / norm(f);
            theta = mod(atan2d(f(1),f(2)), 180);
            amplitude = abs(su(1,1,Np*(k-1)+l,1)) / abs(su(1,1,Np*(k-1)+1,1));
            fprintf('%.2fpx %.2fdeg k=(%.2f,%.2f)px %.2f%% %.2f%%\n', ...
                period, theta, f(1), f(2), norm(f) / cutoff * 100, amplitude*100);
            if verbose_mode >= 1
                figure(4);
                subplot(Na,Np-1,(Np-1)*(k-1)+l-1);
                imshow(fftshift(log(1+abs(S(:,:,l)) .* D)),[]);
                hold on;plot(size(S,1)/2+1+f(2),size(S,2)/2+1+f(1),'yo');hold off;
                title(sprintf('%.2fpx %.2fdeg', period, theta));
                drawnow
            end
        end
        su(:,:,Np*(k-1)+l) = su(:,:,Np*(k-1)+l) .* sw(:,:,Np*(k-1)+l);
    end    
end
% display the shifted components
if verbose_mode >= 1
    figure(5);
    for l=1:L
        subplot(Na,Np,l), imshow(ifftshift(log(1+su(:,:,l))),[]);
        hold on; plot(size(su,1)/2+1,size(su,2)/2+1,'+'); hold off;
    end
    figure(6);
    for k = 1:Na
        idx = Np*(k-1)+1:Np*k;
        if Np == 3
            tmp = fftshift(log(1+su(:,:,idx)));
            tmp = uint8((tmp - min(tmp(:))) ./ (max(tmp(:))-min(tmp(:))) * 255);
            subplot(3,Na,k), imshow(tmp,[0 255])
        end        
        tmp = mean(conj(sw(:,:,idx)) .* su(:,:,idx), 3) ./ (mean(abs(sw(:,:,idx)).^3, 3) + v);
        tmp = real(ifft2(tmp));
        subplot(3,Na,Na+k), imshow(tmp,[])
        subplot(3,Na,2*Na+k), fftshow(tmp)
    end
    figure(7);
    imshow(fftshift(log(1+mean(abs(sw),3))),[])
    title('Sum of shifted OTFs');
end

fprintf('Combine images (%s)\n', method);

su = mean(su, 3) ./ (mean(abs(sw).^3, 3) + v);      
su = real(ifft2(su));


if (zoom == 2)    
    M = fftupscale(fft2(mean(data,3)));    
    dec = real(ifft2(conj(H) .* M./(v + abs(H).^2)));
    avg = real(ifft2(M));
else
    dec = real(ifft2(conj(otf).*fft2(mean(data,3))./(v.^2 + abs(otf).^2)));
    avg = mean(data,3);
end
if padsize > 0     
    su = su(zoom*padsize:zoom*N-zoom*padsize-1,zoom*padsize:zoom*N-zoom*padsize-1);  
    dec = dec(zoom*padsize:zoom*N-zoom*padsize-1,zoom*padsize:zoom*N-zoom*padsize-1);
    avg = avg(zoom*padsize:zoom*N-zoom*padsize-1,zoom*padsize:zoom*N-zoom*padsize-1);
end
end

function I2 = register_int(H1,H2,I1,I2)
% compute a linear regression between the ref I1 and I2 to determine if I2
% should be negated. Why is this needed?
J = (H1 > eps) & (H2 > eps);
if sum(J(:)) > 4
    X = imag(I1(J));
    Y = imag(I2(J));
    p = polyfit(X(:), Y(:), 1);    
    I2 = I2 * sign( p(1));%(I2 - p(2)) / p(1);
    %p(1)
end
end

function f = peak(img,subpixel_localization)
img = fftshift(img);
[~,k] = max(img(:));
[dy, dx] = ind2sub(size(img), k);
if subpixel_localization
    [dx,dy] = fit_quadratic_peak(img,dx,dy,1);
end
%figure(12), clf, imshow(img,[]);hold on; plot(dx,dy,'ro');hold off;pause;
dx = dx - (size(img,2))/2 - 1;
dy = dy - (size(img,1))/2 - 1;
f = [dy;dx];
end

function f = peak2(im1,im2,D,subpixel_localization)
% register im1 and im2 using phase correlation
im2 = im2 - im1;
im1 = fftshift(abs(im1));
im2 = fftshift(abs(im2));
%figure(12);subplot(121),imshow(im1,[]);subplot(122);imshow(im2,[]);pause
r = fftn(im1) .* conj(fftn(im2));
r = ifftshift(D.*real(ifftn((r./abs(r)))));
[~,k] = max(r(:));
[dy, dx] = ind2sub(size(r), k);
if subpixel_localization
    [dx,dy] = fit_quadratic_peak(r,dx,dy,2);
end
%figure(12), clf, imshow(r,[]);hold on; plot(dx,dy,'ro');hold off;pause();
dx = dx - (size(r,2))/2 - 1;
dy = dy - (size(r,1))/2 - 1;
f = [dy;dx];
end

function dest = phasesplit(data)
M = exp(-2*pi*1i*[0 0 0;0 -1 1;0 -2 2]/3);
dest = zeros(size(data));
for r = 1:size(data,1)
    for c = 1:size(data,2)
        v = reshape(data(r, c, :), 3, 1);
        dest(r, c, :) = M \ v;
    end
end
end

function w = weights(n,k,m,s)
[x, y] = meshgrid(-n/2:(n/2-1),-n/2:(n/2-1));
w = sqrt((x-k(2)).^2+(y-k(1)).^2);
w = 1 - m * exp(-w/s);
%w = 1 - m*(sqrt(((x-k(2)).^2+(y-k(1)).^2)/(s^2)) < 1);
w = fftshift(w);
end

function w = distance2(n,x0,y0)
[x, y] = meshgrid(-n/2:(n/2-1),-n/2:(n/2-1));
w = sqrt((x-x0).^2 + (y-y0).^2);
w = fftshift(w);
end

function noise_std = noise_level(data)
tmp = zeros(size(data));
for l = 1:size(data,3)
    tmp(:,:,l) = del2(data(:,:,l));
end
noise_std = 1.48 * mad(tmp(:),1);
end

function [data,Na,Np] = reorderimages(data,otf,verbose_mode)
N = size(data,1);
thetas = zeros(size(data,3),1);
D = otf > 0 & otf < 0.4;
if verbose_mode > 0
    clf;
end
A = mean(data,3);
for l = 1:size(data,3)
    S = log(1+abs(fft2(data(:,:,l) - A)));
    %S = S-imfilter(S,fspecial('gaussian',[21 21], 3),'symmetric');
    S = S .* D;
    f = peak(S,0);
    thetas(l) = mod(atan2d(f(1),f(2)),180);    
end

[~,perm] = sort(round(thetas));

data = data(:,:,perm);
Na = numel(unique(round(thetas)));
Np = size(data,3) / Na;
end

function su = apotome( data )
% img = APOTOME(data)
% reconstruction of apotome images
%
% [1] Neil, M. A. A., Juškaitis, R. & Wilson, T. Method of obtaining optical sectioning by using structured light in a conventional microscope. Optics Letters 22, 1905 (1997).
% [2] O’Holleran, K. & Shaw, M. Polarization effects on contrast in structured illumination microscopy. Optics Letters 37, 4603 (2012).
%
% Jerome Boulanger 2017

avg = mean(data,3);
su = zeros(size(data,1), size(data,2));
  for k = 1:size(data,3) / 3
    idx = 3*(k-1)+1:3*k;    
    a = (data(:,:,idx(1))-data(:,:,idx(2))).^2;
    b = (data(:,:,idx(1))-data(:,:,idx(3))).^2;
    c = (data(:,:,idx(2))-data(:,:,idx(3))).^2;
    %u = sqrt(2) / 3 * sqrt(a+b+c);  
    u = 3/2 * sqrt(a+b+c);  % this give 100% contrast for 
    fprintf(1,'Contrast %d %.2f%%\n', k, 100*mean(u(:)) ./ mean(avg(:)));
    su = su + u;
  end
  su = 3 * su / size(data,3);% + mean(data,3);    
end

