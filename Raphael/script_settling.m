clear all;
close all;

dat= importdata('./data/PSD_DATA.mat');

% itest = 8;
% dat.Sample_id
% d = dat.Diameter{itest}/1000;
% hist_part = dat.Vol_percent{itest};
% 
% vs = 2/9*1.5*9.81*(d/2).^2;


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


figure
plot(vs,hist_part)
Pv  = hist_part/sum(hist_part);
export_fig(gcf,'PVD.png');

c0tot = 0.5;
t = 0:600:7200
ci0 = c0tot * Pv;
z = 0:-0.0001:-0.1

ci = @(ci0,vs,z,t)(ci0*(1-heaviside(z+vs*t)));

%
figure
C = zeros(length(z),length(t));
for k = 1:1:length(t)
    %t(k)/t(end)*100
    Ci = 0;
    for i = 1:length(Pv)
        Ci = ci(ci0(i),vs(i),z,t(k)) + Ci;
    end
     plot(Ci,z,'LineWidth',3)
     hold on
%     xlim([0 2*c0tot]);
%     ylim([min(z) max(z)]);
%     pause(0.1)
   % C(:,k) = Ci';
    
end
set(gcf,'position',[1640 841 560 1217]);
export_fig(gcf,'profiles.png');
