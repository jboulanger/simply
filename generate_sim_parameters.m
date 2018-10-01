function p = generate_sim_parameters(period,Na,Np)
% p = generate_sim_parameters(period,Na,Np)
% 
% Generate a set of parameters for SIM data
%

L = Na*Np;
period = period * ones(Np,Na);
theta = kron((0:Na-1)/Na*180, ones(Np,1)) + 15;
shifts = kron(ones(1,Na), (0:Np-1)'/Np);
amplitudes = 0.1 + 0.05*rand(Np,Na);
for l = 1:L
    p(l).period = period(l);
    p(l).orientation = theta(l);
    p(l).shift = shifts(l);
    p(l).amplitude = amplitudes(l);
end