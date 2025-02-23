clear
% close all

sizeP = 350
exp = {'SiCF500_05.1.8.14','SiCF500_05.2.8.14','SiCF800_05.1.8.14','SiCF800_05.2.8.14','SiCF1000_05.1.8.14','SiCF1000_05.2.8.14','SiC33_05.1.8.16','SiC33_05.2.8.16','SiC50_05.1.8.17','SiC50_05.2.8.17'}
exp = {'CCFZ_05_sat_180'};
vmax = [0.02,0.02,0.01,0.01,0.0025,0.0025,0.02,0.02,0.02,0.02]/35

for dd = 1:1
PVDML = load(['/Users/marcoturrini/Desktop/PVDML_ds_' num2str(dd) '.txt']);
load(['./exp_data/' exp{dd} '.mat']) % load(['/Users/marcoturrini/Desktop/test_set/data_set_' num2str(dd) '.mat']);
t = dat.t;
z = dat.z;
tplt = logspace(log10(t(2)),log10(t(end)),10);
Ptot = PVDML/sum(PVDML);
c0 = mean(dat.C(ceil(0.4*length(z)):ceil(0.6*length(z)),1));
% vref = logspace(0,3,sizeP)/10^3*dat.v(end);
vref = logspace(0,3,sizeP)/10^3*vmax(dd);
ci = @(ci0,vs,z,t)(ci0*(1-heaviside(z+vs*t)));
C_t_z_or = (dat.C-min(min(dat.C)))/c0*0.4;
c0 = 0.5


C = zeros(length(z),length(t));
cnt = 0;
col = lines(100);
figure
bar(log(vref),Ptot,'EdgeAlpha',[0],'FaceAlpha',0.3,'DisplayName',['PVD']);
hold on
% bar(log(dat.v),dat.P,'Facecolor',[0 0.3 1],'EdgeAlpha',[0],'FaceAlpha',0.3,'DisplayName',['PVD ANN ']);

pause(2)
figure
hold on
for k = 1:length(tplt)
    
    [~,itp] = min(abs(tplt(k)-t));
    tp = t(itp);
    cnt = cnt+1;
    %t(k)/t(end)*100
    Ci = 0;
    for i = 1:length(Ptot)
        Ci = ci(c0*Ptot(i),vref(i),z,tp) + Ci;
    end
    C(:,k) = Ci;
    plot(Ci,z,'color',col(cnt,:)); hold on;
    plot(C_t_z_or(:,itp),z,'--','color',col(cnt,:));
end
     xlabel('C [g/L]');
     ylabel('z [m]');
     set(gca,'FontSize',14);

pause(2)
end

set(gcf,'Position',[440   181   415   617])
set(gca,'box','on')