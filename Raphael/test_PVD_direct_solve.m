clear all;
close all;

dat = importdata('dat.mat');

addpath('../Toolboxes/');

t = dat.t;
z = dat.z;
C_t_z = dat.C;

z0 = 0.09;
z1 = 0.089;


vmax = 2e-4;
[vi,Pi] = PVD_direct_solve(t,z,C_t_z,z0,z1,vmax);

figure
bar(vi,Pi);
ylabel('Pi');


c0 = 0.5;
ci = @(ci0,vs,z,t)(ci0*(1-heaviside(z+vs*t)));
stp = 20;
figure
C = zeros(length(z),length(t));
cnt = 0;
col = lines(100);
for k = 1:stp:length(t)
    cnt = cnt+1;
    %t(k)/t(end)*100
    Ci = 0;
    for i = 1:length(Pi)
        Ci = ci(c0*Pi(i),vi(i),z,t(k)) + Ci;
    end
   
    plot(Ci,z,'color',col(cnt,:)); hold on;
    plot(C_t_z(:,k),z,'--','color',col(cnt,:));
end

    



