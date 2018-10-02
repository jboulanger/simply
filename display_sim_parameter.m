function display_sim_parameter(p,pixel_size,wavelength,numerical_aperture)
%
% display_sim_parameter(p,pixel_size,wavelenght,numerical_aperture)
% display_sim_parameter(p)
%
% Display estimated SIM parameters
%
% Jerome Boulanger 2018

fprintf(1,'period\tangle\tampl.(%%)\tshift(%%)\n');
for l = 1:numel(p)
    fprintf(1,'%.3f\t%.2f\t%.3f\t%.2f\n', p(l).period, p(l).orientation,...
        p(l).amplitude*100,p(l).shift*100);
end

if nargin > 1
    r0 = wavelength / (2*numerical_aperture);
    p0 = mean([p(:).period]) * pixel_size;
    fprintf(1,'period %dnm -> %d%% cutoff\n',  round(p0), round(100*r0 / p0));
end

 fprintf(1,'amplitude avg : %.2f%%\n',  mean([p(:).amplitude])*100);
