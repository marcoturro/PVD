t = dat.t;
z = dat.z;
tplt = logspace(log10(t(2)),log10(t(end)),10);
Ptot = PVDMLds1;
c0 = 1;
vref = linspace(0,1,100)*0.05;
ci = @(ci0,vs,z,t)(ci0*(1-heaviside(z+vs*t)));
C_t_z_or = dat.C(3:end-3,:);

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
    plot(C_t_z_or(:,itp),z,'--','color',col(cnt,:));
end
     xlabel('C [g/L]');
     ylabel('z [m]');
     set(gca,'FontSize',14);