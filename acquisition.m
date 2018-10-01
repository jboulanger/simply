%% Reconstruct acquired SIM data
% 
% Jerome Boulanger 2018

clear all

% select a file
[f, d] = uigetfile('../Data/*.tif');
filename = [d, f];
fprintf('Loading ''%s''\n', filename);
metadata = imfinfo(filename);
% load the data in a 3D array
data = zeros(metadata(1).Height,metadata(1).Width,numel(metadata));
for l = 1:numel(metadata)
  data(:,:,l) = double(imread(filename, l));
end

% padding in case images are not square
if size(data,1) < size(data,2)    
    P = (size(data,2)-size(data,1));
    data = padarray(data,[P,0,0],min(data(:)),'post');
end

if size(data,1) > size(data,2)    
    P = (size(data,1)-size(data,2));
    data = padarray(data,[P,0],min(data(:)),'post');
end

% load or define microscopy parameters
pfile = [d strrep(f,'.tif','.mat')];
if exist(pfile,'file')
    load(pfile);
else    
    answer = inputdlg({'Pixel_size','Wavelength','NA'},...
        'Microscope parameter', 1, {'82.9','580','1.49'}); 
    pixel_size = str2double(answer{1});
    wavelength = str2double(answer{2});
    numerical_aperture = str2double(answer{3});
    fprintf('Saving parameters to file ''%s''.\n', pfile);
    save(pfile,'pixel_size','wavelength','numerical_aperture');
end
cutoff = wavelength / (2 * pixel_size * numerical_aperture);
fprintf('px:%.2fnm lambda:%.0fnm NA:%.2f (%.2fpx)\n', ...
    pixel_size, wavelength, numerical_aperture, cutoff);
otf = generate_otf(max(size(data)), cutoff);

%% estimate modulations
tic
phat = estimate_sim_parameters(data,otf);
display_sim_parameter(phat,pixel_size,wavelength,numerical_aperture)
toc
%%  reconstruct the image
tic
zoom = 2;
wiener_parameter = 1;
mask_amplitude = 0.75;
[im,sw] = reconstruct_sim_base(data,phat,otf,zoom, ...
    wiener_parameter,mask_amplitude);
toc

figure(1);
subplot(221), imshow(mean(data,3),[]);
subplot(222), fftshow(mean(data,3),otf);
subplot(223), imshow(im,[0 max(im(:))])
subplot(224), fftshow(im,otf);

%% save the file
ofilename = strrep(filename,'.tif','_sim.tif');
fprintf('Saving file ''%s''\n', ofilename);
imwrite(uint16(max(0,im)),ofilename);
