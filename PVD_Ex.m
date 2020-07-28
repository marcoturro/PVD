clear all;
close all;

in = 'D:\Sediment_exp\Measurements\Marco05.7.7.20\Marco05.7.7.20.mat'
dat = importdata(in);
addpath('./Toolboxes')

t = dat.t;
[~, id_z0] = min(abs(dat.z));
z = dat.z(id_z0:end);
C_t_z = dat.C(id_z0:end,:);

C_t_z = C_t_z - min(min(C_t_z));
 
for i = 1:length(t)
    id_zM = 40;
    C_t_z(1:id_zM,i) = interp1([0 z(id_zM)],[0 C_t_z(id_zM+1,i)],z(1:id_zM));
end



z0 = 0.07;
zz1 = linspace(z0/5,z0*0.99,10);
for lopt = 1:length(zz1)
z1 = zz1(lopt)


vmax = 1e-1;
[vi,Pi] = PVD_direct_solve(t,z,C_t_z,z0,z1,vmax);

figure
bar(log(vi),Pi);
ylabel('Pi');


c0 = 0.5;
ci = @(ci0,vs,z,t)(ci0*(1-heaviside(z+vs*t)));


tplt = logspace(log10(t(2)),log10(t(end)),10);


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

vi_z{lopt} = vi;
Pi_z{lopt} = Pi;


end


vref = vi_z{lopt}
for jj = 1:lopt
    Ps(jj,:) = interp1(vi_z{jj},Pi_z{jj},vref,'pchip');
end

Ptot = mean(Ps,1)/lopt;
Ptot = Ptot/sum(Ptot);
figure
bar(log(vref),Ptot);
ylabel('Pi');

c0 = 0.5;

tplt = logspace(log10(t(2)),log10(t(end)),10);

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
    for i = 1:length(Ptot)
        Ci = ci(c0*Ptot(i),vref(i),z,tp) + Ci;
    end
   
    plot(Ci,z,'color',col(cnt,:)); hold on;
    plot(C_t_z(:,itp),z,'--','color',col(cnt,:));
end
