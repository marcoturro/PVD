clear all;
close all;
tic
Disp = 0;
addpath('./Toolboxes')

in = 'E:\Raphael\MIT\sediment_model\GIT\PVD\Raphael\created_data\data_set_33.mat';
data = importdata(in);

ci = @(ci0,vs,z,t)(ci0*(1-heaviside(z+vs*t)));



try
    t = data.t;
    z = data.z; % attention untis of z
    C_t_z = data.C;
    [m, ~] = size(C_t_z);
    if m ~= length(z)
        C_t_z = C_t_z';
    end
    if isfield(data.dat,'P')
        Pr = data.P;
        vr = data.v;
        %vr = vr./(maxz/maxt);
        realP = 1;
    else
        realP = 0;
    end
catch
    
end


zz1 = [0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,0.95,0.98];
for lmnop = 1:length(zz1)
z1 = zz1(lmnop);
z0 = 1;


z1dim = z1 * max(abs(z));
z0dim = z0 * max(abs(z));
%sprintf(' z1 = %3.3f m ;  z0 = %3.3f m',z1dim,z0dim)

c0 = mean(C_t_z(:,1));
cN = mean(C_t_z(:,end));

%C_v = C_t_z./c0.*(1-C_t_z./c0);

[~, z1id] = min( abs( z + z1dim) );
[~, z0id] = min( abs( z + z0dim) );

Ci = zeros(1,length(t)) ; %Ciwvl = Ci;

imax = min(find(1-abs(C_t_z(z1id,:)/c0) > 0.01))-1; %is this the best way to do it?

imin = max(find(abs(C_t_z(z0id,:)/c0) > 0.01));
if isnan(imin)
    imin = length(t);
end
for i = 1:length(t)
    Ci(i)     = trapz(z(z1id:z0id),-C_t_z(z1id:z0id,i));
    
end
tvmin = t(imin);
Pvmin = abs(Ci(imin)/c0);
tvmax = t(imax);
dv = 5  * z1*(1/t(imax+1)-1/t(imax));
Pvmax = 1-abs(Ci(imax)/c0);


Ci = Ci./max(Ci);

dat.t =t;
dat.z = z;
dat.Ci = Ci;
dat.c0 = c0;
dat.z0 = z0;
dat.z1 = z1;
dat.C_t_z = C_t_z;
dat.tvmin = tvmin;
dat.tvmax = tvmax;
dat.dv = dv;
dat.Pvmin =Pvmin;
dat.Pvmax = Pvmax;
dat.Disp = Disp;
%= %struct('t',t,'z',z,'Ci',Ci,'c0',c0,'z0',z0dim,'z1'...
% ,z1dim,'Disp',Disp,'C_t_z',C_t_z);%,'C_v',C_v); % Chose Ci or Ciwvl

[vs, P, ERR] = PVD_solve(dat);
length(P)




P = P/sum(P);

Nv = length(P);

stp = fix(length(t)/20);

Cmat = produce_data(P',vs,t(1:stp:end),z)';

Ctest = C_t_z(:,1:stp:end);
err(lmnop) = 1/length(Cmat(:))*sum((Cmat(:)-Ctest(:)).^2);

plt = 0;
if plt
figure(lmnop)
plot(Cmat,z,'+'); hold on;
hold on
title('Concentration profiles')
xlim([0 c0*1.1]); ylim([min(z),max(z)])
plot(C_t_z(:,1:stp:end),z);
%scatter([C_t_z(boundD,k) C_t_z(boundU,k)],[z(boundD) z(boundU)])
xlabel('C [g/L]'); ylabel('z [cm]')
legend('reconstructed data','input data');
set(gca,'FontSize',18)

if realP
    figure
    plot(vs,P); hold on;
    plot(vr,Pr);
end
end

end

figure
plot(zz1,err);
toc