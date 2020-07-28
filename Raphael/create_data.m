clear all;
close all;

Nv = 4000;
Ns = 20
sig = logspace(-6,-3.9,Ns);
var = 2e-6*ones(size(sig));
var(8:end) = 1e-5;
var(10) = 1e-6;
Np = floor(logspace(6,6.2,Ns));

y = 0
for i=1:length(sig)
yv = randraw('normal', [sig(i),var(i)], Np(i)); %Creates a random gaussian distribution of particles centered on diameter dp50 with variance var_dp
y = [y yv'];
end

[hist_part,vs] = hist(y,Nv); %Calculates the histogram of particle size with Nd bins, and returns the histogram and the mean particle size for each bin
%hist_part = fliplr(hist_part)
%hist_part(end) = 5e4
% vs = logspace(-6,-3,Nv);
% hist_part = ones(size(vs))


Pv  = hist_part/sum(hist_part);

c0tot = 0.5;
t = logspace(-5,5,1000);
ci0 = c0tot * Pv;
z = 0:-0.0005:-0.1;

ci = @(ci0,vs,z,t)(ci0*(1-heaviside(z+vs*t)));


C = zeros(length(z),length(t));
for k = 1:1:length(t)
    %t(k)/t(end)*100
    Ci = 0;
    for i = 1:length(Pv)
        Ci = ci(ci0(i),vs(i),z,t(k)) + Ci;
    end
    C(:,k) = Ci;
end

dat.C = C;
dat.z = z;
dat.t = t;

save('dat.mat','dat');