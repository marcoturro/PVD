clear; close all; tic
addpath('./Toolboxes')

pltC = 1; %display the concentration profiles
pltP = 1; % display the PVDs
pltE = 1; % display the 
No = 1000;

in = './exp_data/photo2.7.13.mat';
% in = './Raphael/created_data/data_set_41.mat';
dat = importdata(in);

try
    vmax = dat.v(end);
    [C_t_z, z, t, vr, Pr] = data_norm(dat,vmax);
catch
    vmax = input('vmax? \n'); %there is two vmax
    [C_t_z, z, t, vr, Pr] = data_norm(dat,vmax);
    realP = 0;
end


zz1 = [0.3,0.4,0.5,0.6,0.7,0.8,0.9,0.95,0.96,0.97,0.98,0.99];
%zz1 = 0.3
%c0 = 1; 
c0 = mean(C_t_z(:,1));

for lmnop = 1:length(zz1)
z1 = zz1(lmnop);
z0 = 1;
[~, z1id] = min( abs( z + z1) );
[~, z0id] = min( abs( z + z0) );

Ci = zeros(1,length(t)) ; 

% this section aims to find the spatial limits within which the relevant
% information is found.

for i = 1:length(t)
    Ci(i)     = trapz(z(z1id:z0id),-C_t_z(z1id:z0id,i));
end
% 
Ci = wdenoise(Ci,8);
Ci = movmean(Ci,3);
order = 3; frame = 15;
Ci = sgolayfilt(Ci,order,frame);

imax = find(abs(diff(Ci)./diff(t)) > 1e-8,1) - 1;
imin = find(abs(diff(Ci)./diff(t)) > 1e-8, 1, 'last' );

if isnan(imin)
    imin = length(t);
end

if imax < 5 
    imax = 10;
end

tvmin = t(imin);
Pvmin = abs(Ci(imin)/c0);
tvmax = t(imax);
Pvmax = abs(Ci(imax)/c0);

solve.t =t;
solve.z = z;
solve.Ci = Ci;
solve.c0 = c0;
solve.z0 = z0;
solve.z1 = z1;
solve.tvmin = tvmin;
solve.tvmax = tvmax;
solve.Pvmin = Pvmin;
solve.Pvmax = Pvmax; 


[vs, P] = PVDS(solve);
P%  = P.*[vs(1), diff(vs)];
length(P)

P = P/max(P);
P = P/sum(P);
Nv = length(P);
stp = fix(length(t)/7);
Cmat = produce_data(P',vs,t(1:stp:end),z)';
Ctest = C_t_z(:,1:stp:end);
err(lmnop) = 1/length(Cmat(:))*sum((Cmat(:)-Ctest(:)).^2);


if pltC
figure
hold on
for col = 1:length(Cmat(1,:))
R = rand; G = rand; B = rand;
plot(Cmat(:,col),z,'-+','color',[R,G,B]); 
plot(Ctest(:,col),z,'--','color',[R,G,B]);
end
title(['Concentration profiles z1 = ' num2str(z1)])
xlim([0 1.1]); ylim([min(z),max(z)])
xlabel('C [g/L]'); ylabel('z')

% legend('reconstructed data','input data');
set(gca,'FontSize',18)

if pltP
    figure
    plot(vs,P); hold on;
    plot(vs,movmean(P,10));
    plot(vr,Pr);
end
end

end

if pltE
    figure
    plot(zz1,err);
end

toc