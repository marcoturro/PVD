clear all;
close all;

dat = importdata('photo1.7.13.mat');

t = dat.t;
z = dat.z;
C_t_z = dat.C;
zact = -0.02;
[~,izact] = min(abs(zact-z));
z = z(izact:end)-zact;
C_t_z = C_t_z(izact:end,:);


z0 = 0.05;
z1 = 0.045;


vmax = 1e-1;
[vi,Pi] = PVD_direct_solve(t,z,C_t_z,z0,z1,vmax);

figure
bar(log(vi),Pi);
ylabel('Pi');


c0 = 0.5;
ci = @(ci0,vs,z,t)(ci0*(1-heaviside(z+vs*t)));


tplt = logspace(log10(t(2)),log10(t(end)),20);


figure
C = zeros(length(z),length(t));
cnt = 0;
col = lines(100);
for k = 1:length(tplt)
    
    [~,itp] = min(abs(tplt(k)-t));
    tp = t(itp);
    cnt = cnt+1;
    %t(k)/t(end)*100
    Ci = 0;
    for i = 1:length(Pi)
        Ci = ci(c0*Pi(i),vi(i),z,tp) + Ci;
    end
   
    plot(Ci,z,'color',col(cnt,:)); hold on;
    plot(C_t_z(:,itp),z,'--','color',col(cnt,:));
end

    


