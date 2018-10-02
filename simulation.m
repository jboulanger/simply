% Simulation of structured illumination microscopy imaging and
% reconstruction.
%
% Jerome Boulanger 2018

clear all
filename = 'periodic_pattern.tif';
period = 4.75; % modulation period (before zoom)
cutoff = 4; % cutoff frequency
zoom = 2; % zoom factor (1 or 2)
noise_std = 2; % Gaussin noise level
wiener_parameter = 0.03; % wiener parameter
mask_amplitude = 0.75; % mask amplitude
obj = double(imread(filename)); % load the test image
p = generate_sim_parameters(period,3,3); % generate SIM parameters
m = generate_sim_modulation(size(obj), p); % generate a set of modulation
otf = generate_otf(size(obj,1), cutoff); % compute an OTF
[data,otf] = generate_sim_data(obj,m,otf,zoom,noise_std); % generate SIM data
tic;
phat = estimate_sim_parameters(data,otf); % estimate SIM parameters
display_sim_parameter(phat); % display a summary of the parameters
mhat = generate_sim_modulation(size(obj), phat); % and compute the modulation
[im,sw] = reconstruct_sim_base(data,phat,otf,zoom,wiener_parameter,mask_amplitude);
toc;
% compute the root mean square error
c = polyfit(obj(:),im(:),1);
im = (im - c(2)) / c(1);
avg = mean(data,3);
avg = real(ifft2(ifftshift(padarray(fftshift(fft2(avg)),size(avg)/2))));
c = polyfit(obj(:),avg(:),1);
avg = (avg - c(2)) / c(1);
fprintf(1,'Error AVG RMS:%f\n', sqrt(mean((avg(:)-obj(:)).^2))/mean(obj(:)));
fprintf(1,'Error SIM RMS:%f\n', sqrt(mean((im(:)-obj(:)).^2))/mean(obj(:)));
% display results
figure(1), clf;
subplot(221), imshow(avg,[]);
subplot(223), fftshow(avg,otf);
subplot(222), imshow(im,[]);
subplot(224), fftshow(im,otf);
