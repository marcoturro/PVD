clear
% close all;

in = ['./exp_data/SiCF1000_05.2.8.14.mat']
% in = './Marco/created_data/data_set_1.mat'
dat = importdata(in);
dat.z = dat.z;
addpath('./Toolboxes')
plt = 1; 
mPvd = 1;
vmax = 0.5e-1;

t = dat.t;
[~, id_z0] = min(abs(dat.z));
z = dat.z(id_z0:end);

if length(dat.C(1,:)) ~= length(t)
    Cz = dat.C';
    C_t_z = Cz(id_z0:end,:);
else
    C_t_z = dat.C(id_z0:end,:);
end

C_t_z = C_t_z(1:end,:) - min(min(C_t_z));
C_t_z_or = C_t_z;
C_t_z = wdenoise(C_t_z,4);

for i = 1:length(t)
    id_zM = 20;
    C_t_z(1:id_zM,i) = interp1([0 z(id_zM)],[0 C_t_z(id_zM+1,i)],z(1:id_zM));
end

c0 = mean(C_t_z(200:end-200,1));
ci = @(ci0,vs,z,t)(ci0*(1-heaviside(z+vs*t)));

zz0 = 0.02:0.001:0.09;
zz0 = 0.06;

clear Ps vs

for lopt = 1:length(zz0)
z0 = zz0(lopt);
z1 = z0*0.99;

[vi,Pi,pi] = PVD_direct_solve(t,z,C_t_z,z0,z1,vmax);
    


    tplt = logspace(log10(t(2)),log10(t(end)),10);


%     Prof_fig = figure
    C = zeros(length(z),length(t));
    cnt = 0;
    col = lines(100);
%     for k = 1:length(tplt)
% 
%         [~,itp] = min(abs(tplt(k)-t));
%         tp = t(itp);
%         cnt = cnt+1;
%         Ci = 0;
%         for i = 1:length(Pi)
%             Ci = ci(c0*Pi(i),vi(i),z,tp) + Ci;
%         end
% 
%         plot(Ci,z,'color',[col(cnt,:) 0.1]); hold on;
%         plot(C_t_z_or(:,itp),z,'--','color',col(cnt,:));
%     end
%     xlabel('C [g/L]');
%     ylabel('z [m]');
%     set(gca,'FontSize',14);


vi_z{lopt} = vi;
Pi_z{lopt} = pi;


end

if mPvd == 1
vref = vi_z{lopt};

for jj = 1:lopt
    Ps(jj,:) = interp1(vi_z{jj},Pi_z{jj},vref,'pchip');
end

Ptot = mean(Ps,1);
Ptot = Ptot.*[vref(1) diff(vref)];
Ptot = Ptot/sum(Ptot);

PVD_fig = figure;
bar(log(vref),Ptot,'Facecolor',[0 0.3 1],'FaceAlpha',0.3,'DisplayName',['z0 = ' num2str(z0)]);
ylabel('Pi');
xlabel('log(v) [m/s]');
set(gca,'FontSize',14);


tplt = logspace(log10(t(2)),log10(t(end)),10);

% figure(Prof_fig)
C = zeros(length(z),length(t));
cnt = 0;
col = lines(100);

% for k = 1:length(tplt)
%     
%     [~,itp] = min(abs(tplt(k)-t));
%     tp = t(itp);
%     cnt = cnt+1;
%     %t(k)/t(end)*100
%     Ci = 0;
%     for i = 1:length(Ptot)
%         Ci = ci(c0*Ptot(i),vref(i),z,tp) + Ci;
%     end
%     
%     plot(Ci,z,'color',col(cnt,:)); hold on;
%     plot(C_t_z_or(:,itp),z,'--','color',col(cnt,:));
% end

end