function [im,sw,C] = reconstruct_sim_base(data,p,otf,zoom,wiener_parameter,mask_amplitude)
%  im = reconstruct_sim_base(data,p,otf,zoom,wiener_parameter)
%
%  Structured illumination image reconstruction using Gustafsson
%  reconstruction method.

% Input :
% - data is a 3d arrays with angles x phases in the 3d dimension
% - p is a array of structures period [px], orientation [deg], 
%   phase [fraction of 2*pi], amplitude
% - otf is optional and used to normalized amplitudes (centered in (1,1))
% - zoom can be 1 or 2
% - wiener_parameter is the regularization parameters used for the wiener
%   filter
%
% Output:
% - im is the reconstructed SIM image
% - sw is the fourier weights used in the reconstruction
% - C are the individual (fft) components after demodulation
%
% Jerome Boulanger 2018

% find the discrete set of angle
tol = 5;
theta = unique(round([p(:).orientation]/tol)*tol);

[x,y] = meshgrid(0:zoom*size(data,2)-1,0:zoom*size(data,1)-1);

% pad the otf if zooming
if zoom == 2
    w0 = ifftshift(padarray(fftshift(otf),[size(otf,1)/2,size(otf,2)/2],'both'));
else
    w0 = otf;
end

mask = 1 - mask_amplitude * w0;

swu = zeros(zoom * size(data,1), zoom * size(data,2));
sw = zeros(zoom * size(data,1), zoom * size(data,2));

if nargout > 2
    C = zeros(size(w0,1),size(w0,2),size(data,3));
end

fprintf('    Period \t|   Orientation\t|  Amplitude\n')
fprintf('mean\t std\t|  mean\t std \t| mean\t std\n')

for k = 1:numel(theta)    
    
    % compute the indices corresponding to 1 angle
    idx = find(abs([p(:).orientation]-theta(k)) < tol);   
    thetak = mean([p(idx).orientation]);
    
    if numel(idx) == 3
    % compute average of period and orientation
    
    pk = zoom * mean([p(idx).period]);    
    fprintf('% 2.2f\t% 2.2f\t| %-3.1f\t% 2.2f\t| %.2f\t%.2f\n', ...
        pk,     log10(std([p(idx).period])), ...
        thetak, log10(std([p(idx).orientation])), ...
        mean([p(idx).amplitude]), log(std([p(idx).amplitude])));
    
    z = 2 * pi * (cosd(thetak).*x + sind(thetak).*y) / pk;
    
    % in fourier space
    kx = round(size(sw,2) * cosd(thetak) / pk);
    ky = round(size(sw,1) * sind(thetak) / pk);
    
    % build the phase matrix and inverse it
    M = zeros(numel(idx),numel(idx));    
    for j = 1:numel(idx);
        phi = 360 * p(idx(j)).shift;
        a = p(idx(j)).amplitude;        
        M(j,:) =  0.5 * [1, a * cosd(phi), -a * sind(phi)];  
    end
    M = inv(M);
    
    % phases splitting and demodulation
    for i = 1:numel(idx)
        u = zeros(size(data,1), size(data,2));
        w = zeros(size(w0));
        % split components
        for j = 1:numel(idx)
            u = u + M(i,j) * data(:,:,idx(j));
            w = w + w0;
        end
               
        % upscale the image and mask 
        if mask_amplitude ~= 0 || zoom ~= 1
            u = fft2(u);
            if zoom == 2
                u = fftshift(u);
                u = padarray(u, [size(u,1)/2,size(u,2)/2], 'both');
                u = ifftshift(u);
            end
            %if i == 2 || i == 3
                u = mask .* u;
                w = mask .* w;
            %end
            u = ifft2(u);
        end
        
        % demodulation and fourier transform 
        switch i            
            case 2                
                w = 0.5 * (circshift(w,[ky,kx]) + circshift(w,[-ky,-kx]));
                u = u .* cos(z);
            case 3
                w = 0.5 * (circshift(w,[ky,kx]) + circshift(w,[-ky,-kx]));                
                u = u .* sin(z);                  
        end          
        
        u = conj(w) .* fft2(u);
        w = conj(w) .* w;
        
        if nargout > 2            
            C(:,:,i+numel(idx)*(k-1)) = u;
        end                
        swu = swu + u;
        sw = sw + w;    
    end
    else 
        fprintf('No enough phases for orientation %.2f\n', thetak);
    end
end

% wiener filter reconstruction
im = real(ifft2(swu ./ (sw + wiener_parameter^2)));


