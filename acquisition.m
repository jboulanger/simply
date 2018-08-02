%% Reconstruct acquired SIM data

clear all
warning off
% select a file
[f, d] = uigetfile('../Data/*.tif');
%%
filename = [d, f];
fprintf('Loading ''%s''\n', filename);
for l = 1:numel(imfinfo(filename))
  data(:,:,l) = double(imread(filename, l));
end
data = data - median(data(:));
%%
if (strfind(filename,'CFK') > 0)
    pixel_size = 83.3;
    wavelength = 610;
    numerical_aperture = 1.2;
elseif (strfind(filename,'Curie') > 0)
    pixel_size = 65;
    wavelength = 500;
    numerical_aperture = 1.49;
elseif (strfind(filename,'SLM-SIM') > 0)
    pixel_size = 90;
    wavelength = 690;
    numerical_aperture = 1.2;  
elseif (strfind(filename,'Zeiss') > 0)
    pixel_size = 65;
    wavelength = 520;
    numerical_aperture = 1.30;    
else
    pixel_size = 82.9;
    wavelength = 580;
    numerical_aperture = 1.49;
end

cutoff = wavelength / (2 * pixel_size * numerical_aperture);
fprintf('px:%.2fnm lambda:%.0fnm NA:%.2f > %.2fpx\n', pixel_size, wavelength, numerical_aperture, cutoff);
otf = gen_otf(size(data,1), cutoff);
zoom = 2;
wiener_factor = 0.05;
padding_level = 0;
prefilter = false;
notchfilter = 1;
subpixel_localization = true;
verbose_mode = 0;

% reconstruct
tic
[img,dec,avg] = simrec(data,otf,zoom,wiener_factor,padding_level,prefilter,notchfilter,subpixel_localization,verbose_mode);
toc

% display
figure(1);
otf = fftupscale(otf);
subplot(231), imshow(sqrt(avg),[]);title('WF');
subplot(234), fftshow(avg,otf);title('|FFT|');
subplot(232), imshow(sqrt(dec),[]);title('Deconvolution');
subplot(235), fftshow(dec,otf);title('|FFT|');
subplot(233), imshow(sqrt(img),[]);title('SIM');
subplot(236), fftshow(img,otf);title('|FFT|');
figure(3), clf, imshow(max(0,dec),[])
figure(2);clf, imshow(max(0,img),[])

% save 
img = img / max(img(:)) * (2^16-1);
dec = dec / max(dec(:)) * (2^16-1);
avg = avg / max(avg(:)) * (2^16-1);
imwrite(uint16(img),strrep(filename,'.tif','_sim.tif'));
imwrite(uint16(dec),strrep(filename,'.tif','_dec.tif'));
imwrite(uint16(avg),strrep(filename,'.tif','_avg.tif'));