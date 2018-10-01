function p = estimate_sim_parameters(data,otf,visu)
%
% p = estimate_sim_param(data,otf,visu)
% p = estimate_sim_param(data,otf)
% p = estimate_sim_param(data)
%
% Estimate modulation in structured illumination microscopy data
%
% Input :
% - data is a 3d arrays with angles x phases in the 3d dimension
% - otf is optional and used to normalized amplitudes (centered in (1,1))
% Output :
% - p is a array of structures period [px], orientation [deg], phase [fraction of 2*pi], amplitude
%
%
% The modulation period should be between 2 and 6 pixels
% 
% Jerome Boulanger & James Manton 2018

if nargin < 3
    visu = false;
end

%
% estimate the period and orientation
if visu
    figure(1);clf;
end

avg = mean(data,3);
[x,y] = meshgrid(0:3*size(data,2)-1,0:3*size(data,1)-1);
D = sqrt((x - size(x,1)/2).^2 + (y - size(y,1)/2).^2);
D = double(D > 0.15 * size(D,1) ) .* double(D < 0.5 * size(D,1));
scales = logspace(2,8,20); % dft oversampling factors
for l = 1:size(data,3)    
    I = tape(del2((data(:,:,l)-avg).*avg),  0.75);    
    r0 = ifft2( I  , 3 * size(I,1), 3 * size(I,2) );        
    r0 = abs(fftshift(r0));    
    r0 = r0 .* D;    
    [ky, kx] = find(r0 == max(r0(:)),1);
    [kx, ky] = fit_quadratic_peak(r0,kx,ky,1);    
    kx = (kx - size(r0,2) / 2 - 1) * size(data,2) / size(r0,2);
    ky = (ky - size(r0,1) / 2 - 1) * size(data,1) / size(r0,1);    
    X = [];
    Y = [];
    for j = 1:numel(scales)
        r1 = abs(dft(I,[ky;kx], scales(j), 255));
        [dy, dx] = find(r1 == max(r1(:)),1);     
        [dx, dy] = fit_quadratic_peak(r1,dx,dy,10);   
        kx = kx - (dx - size(r1,2) / 2 - 1) / scales(j);
        ky = ky - (dy - size(r1,1) / 2 - 1) / scales(j); 
        if visu
            X(j) = kx;
            Y(j) = ky;
            figure(1);
            subplot(221);
            imshow(r0,[]); 
            hold on; 
            plot(3*kx+size(r0,2)/2+1,3*ky+size(r0,1)/2+1,'ro');
            hold off;
            subplot(222);
            imshow(r1,[]); 
            hold on; 
            plot(dx,dy,'ro');
            plot(size(r1,2)/2+1,size(r1,1)/2+1,'g+');
            plot(X,Y);
            hold off;
            drawnow;
            subplot(223);
            plot(X,Y); 
            axis ij;
            xlabel('X');
            ylabel('Y');
            title(sprintf('iteration %d', j))
        end
    end
    p(l).period = size(data,1) / norm([kx;ky]);
    p(l).orientation = mod(atan2d(ky,kx),180);    
end

% estimate phase and amplitude
[x,y] = meshgrid(0:size(data,2)-1,0:size(data,1)-1);
% define  weigths
w = sqrt(abs(reshape(mean(data,3),[size(data,1)*size(data,2),1])));
for l = 1:size(data,3)
    kx = size(data,2) * cosd(p(l).orientation) / p(l).period;
    ky = size(data,1) * sind(p(l).orientation) / p(l).period;
    z = 2 * pi * (kx .* x + ky .* y) / size(data,1);    
    b = reshape(data(:,:,l),[size(data,1)*size(data,2),1]);
    if nargin < 2
        a = 1;
    else
        a = otf(round(norm([kx;ky])));        
        if a == 0
            a = 1;
        end
    end        
    A = [a.*w.*cos(z(:)) -a.*w.*sin(z(:)) w.*ones(numel(z),1)];
    v = (A'*A)\A'*b;       
    p(l).shift = mod(atan2d(v(2),v(1)),360)/360;    
    p(l).amplitude = norm(v(1:2))./v(3);
end




