clear all;
close all;

dat = importdata('photo2.7.13.mat');
addpath('../Toolboxes')
t = dat.t;
z = dat.z;
C_t_z = dat.C;
zact = -0.026;
[~,izact] = min(abs(zact-z));
z = z(izact:end)-zact;

C_t_z = C_t_z(izact:end,:);

a = 1.05;




z0 = 0.05;
z1 = z0/a;

vmax = 1e-2;


%%%

a = z0/z1;
vmin = z0/t(end);
vi(1) = vmin;
i=1;
while vi(i)<vmax
    vi(i+1) = vi(i)*a;
    i = i+1;
end


%%%
options = optimset('Display','iter','TolX',1e-9,'TolFun',1e-9);
pars0 = [1e-7,1];
kappa = @(v,pars)(pars(1)./v+pars(2));
para = fminsearch(@(pars)get_err(t,z,C_t_z,z0,z1,vmax,kappa,pars),pars0,options);

[vi,Pi] = PVD_direct_solve(t,z,C_t_z,z0,z1,vmax,kappa,para);

NZ = length(z);

%
% figure
% bar(log(vi),Pi);
% ylabel('Pi');


c0 = 0.5;
ci = @(ci0,vs,z,t)(ci0*(1-heaviside(z+vs*t)));
tplt = logspace(log10(t(2)),log10(t(end)),20);


C = zeros(length(z),length(t));
cnt = 0;
col = lines(100);

ierr1 = floor(NZ/4);
ierr2 = floor(3*NZ/4);

figure
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


function ERR= get_err(t,z,C_t_z,z0,z1,vmax,kappa,pars)

[vi,Pi] = PVD_direct_solve(t,z,C_t_z,z0,z1,vmax,kappa,pars);

NZ = length(z);

%
% figure
% bar(log(vi),Pi);
% ylabel('Pi');


c0 = 0.48;
ci = @(ci0,vs,z,t)(ci0*(1-heaviside(z+vs*t)));
tplt = logspace(log10(t(3)),log10(t(end)),10);


C = zeros(length(z),length(t));
cnt = 0;
col = lines(100);

ierr1 = floor(NZ/4);
ierr2 = floor(3*NZ/4);

%figure
for k = 1:length(tplt)
    
    [~,itp] = min(abs(tplt(k)-t));
    tp = t(itp);
    cnt = cnt+1;
    %t(k)/t(end)*100
    Ci = 0;
    for i = 1:length(Pi)
        Ci = ci(c0*Pi(i),vi(i),z,tp) + Ci;
    end
    
    err(k) = sqrt(sum(C_t_z(ierr1:ierr2,itp)-Ci(ierr1:ierr2)').^2);
    % plot(Ci,z,'color',col(cnt,:)); hold on;
    % plot(C_t_z(:,itp),z,'--','color',col(cnt,:));
end

ERR = sum(err);

end








