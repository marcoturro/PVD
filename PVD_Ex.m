clear
close all;

in = './exp_data/Marco05.7.7.20.mat'
dat = importdata(in);
addpath('./Toolboxes')
plt = 0; 

t = dat.t;
t = t+10
[~, id_z0] = min(abs(dat.z));
z = dat.z(id_z0:end);
C_t_z = dat.C(id_z0:end,:);

C_t_z = C_t_z - min(min(C_t_z));
 
for i = 1:length(t)
    id_zM = 10;
    C_t_z(1:id_zM,i) = interp1([0 z(id_zM)],[0 C_t_z(id_zM+1,i)],z(1:id_zM));
end

c0 = mean(C_t_z(:,1));
ci = @(ci0,vs,z,t)(ci0*(1-heaviside(z+vs*t)));

for z0 = 0.005:0.01:0.1
zz1 = linspace(z0/100,z0*0.98,10);
clear Ps vs
for lopt = 1:length(zz1)
z1 = zz1(lopt)


vmax = 1e-1;
[vi,Pi,pi] = PVD_direct_solve(t,z,C_t_z,z0,z1,vmax);

if plt == 1
    figure
    bar(log(vi),Pi);
    ylabel('Pi');


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
end


vi_z{lopt} = vi;
Pi_z{lopt} = pi;


end



vref = vi_z{lopt}
for jj = 1:lopt
    Ps(jj,:) = interp1(vi_z{jj},Pi_z{jj},vref,'pchip');
end

Ptot = mean(Ps,1)/lopt;
Ptot = Ptot.*[vref(1) diff(vref)];
Ptot = Ptot/sum(Ptot)

figure(2)
hold on
bar(log(vref),Ptot,'FaceAlpha',0.3,'DisplayName',['z0 = ' num2str(z0)]);
hold off
ylabel('Pi');

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
end
legend show
