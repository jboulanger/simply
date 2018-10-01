function m = generate_sim_modulation(dims, p, type, rounded)
% m = generate_sim_modulation(dims, p, type, rounded)
% m = generate_sim_modulation(dims, p, type)
% m = generate_sim_modulation(dims, p)
%
% Generate structured illumination modulation
%
% Jerome Boulanger 2018
if numel(dims) == 1
    dims = [dims dims];
end
if nargin < 3
    type = 1;
end
if nargin < 3
    rounded = false;
end

L = numel(p);
[x,y] = meshgrid(0:dims(2)-1,0:dims(1)-1);
m = zeros(dims(1),dims(2),L);

for l = 1:L
    kx = cosd(p(l).orientation) / p(l).period;
    ky = sind(p(l).orientation) / p(l).period;
    if rounded == true
        kx = round(dims(2) * kx) / dims(2);
        ky = round(dims(1) * ky) / dims(1);
    end
    m(:,:,l) = kx*x + ky*y + p(l).shift;
    if type==0
        m(:,:,l) = double(mod(m(:,:,l), 1) < 2/3);
    else
        m(:,:,l) = 0.5*(1+p(l).amplitude * cos(2*pi*m(:,:,l)));
    end
end