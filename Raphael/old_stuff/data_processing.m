clear all;
close all;

nstp = 10000;

tfle = importdata('created_data/time.mat');
z = importdata('created_data/depth.mat');


z0 = 0.4;
z1 = 0.1;
a = z1/z0;

[~,iz0] = min((abs(z+z0)));
[~,iz1] = min((abs(z+z1)));

calc = 0;
if calc
    for i = 1:nstp
        i
        dat = importdata(['created_data/Cdata_' num2str(i,'%05d') '.mat']);
        ct(i) = -trapz(z(iz1:iz0),dat(iz1:iz0));
    end
    dat_ct.ct = ct;
    dat_ct.t = tfle(1:nstp);
    dat_ct.z0 = z0;
    dat_ct.z1 = z1;
    save('created_data/ct.mat','dat_ct');
else
    
    dat_ct = importdata('created_data/ct.mat');
    ct = dat_ct.ct;
end


figure
plot(tfle(1:nstp),ct);

vmin = 0.4*10^(-3);%z0/tfle(nstp);
%ddv = z0/(tfle(2)-tfle(1));
vmax  = 2*10^(-3);

ts = z0/vmax;

dt = tfle(2)-tfle(1);

tnxt = ts+dt;

dv = 20*z0*(1/ts-1/tnxt);


vfk = vmin:dv:vmax;%vmin:0.01:vmax;
vfk = fliplr(vfk);

t = z0./vfk;

cint = interp1(tfle(1:nstp),ct,t);



ddcdtt = my_2FD_non_uniform(t,cint);
inan = find(isnan(ddcdtt));
ddcdtt(inan) = 0;
%ddcdtt_s = my_d2dx(t,cint_s);

figure
plot(t,cint);

figure
plot(ddcdtt); hold on;

v = z0./t;

v_star = v*a;
A = zeros(length(t));

for i = 1:length(v)
    idx = min(find(v<v_star(i)));
    
    coef1 = abs(v_star(i) - v(idx));
    coef2 = abs(v_star(i) - v(idx-1));
    dv = abs(v(idx)- v(idx-1));
    if(isempty(dv))
        coef1 = 0; %coef1/dv;
        coef2 = 0;%coef2/dv;
    else
        coef1 =  coef1/dv;
        coef2 = coef2/dv;
    end
    
    A(i,idx) = - coef2 * a^2;
    A(i,idx-1) = - coef1 * a^2;
    A(i,i) = 1;
    
end
A(1,2) = 0;


rhs = (ddcdtt./v.^3)';
rhs(1) = 0;
rhs(end) = 0;
P = linsolve(A,rhs);
figure
%plot (1000*vs,hist_part/max(hist_part)*max(P)); hold on;
plot(1000*v,P,'+');
ylabel('$P(v)$');
xlabel('$v\ (mm/s)$')
legend('original','reconstructed');
%export_fig(gcf,'P_v_reconstructed.png','-m2');


c0tot = 0.5;
P(1) = 0;
P = c0tot*P / sum(P);

figure
plot(v,P); hold on;

% 
% col= lines(100);
% cnt = 0;
% for k = 1:10:length(t)
%     cnt = cnt+1;
%     C_rebuilt = 0;
%     for i = 1:length(P)
%         C_rebuilt = ci(P(i),v(i),z,t(k)) + C_rebuilt;
%     end
%     
%     C = 0;
%     for i = 1:Nv
%         C = ci(ci0(i),vs(i),z,t(k)) + C;
%     end
%     
%     figure(10)
%     plot(C_rebuilt,z,'-','color',col(cnt,:)); hold on;
%     plot(C,z,'--','color',col(cnt,:)); hold on;
%     
% end
% 
% 
% 
% 
% 
% 
% 
