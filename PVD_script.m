clear; close all; tic
addpath('./Toolboxes')

pltC = 1; %display the concentration profiles
pltP = 0; % display the PVDs
pltE = 0; % display the 
No = 1500;

in = 'C:\Users\ENDLab_W10\DSLR\05c\05c.mat';
dat = importdata(in);

% step 1 normalise the data 

dat.C = dat.C-min(dat.C(:,end));
try
    vmax = dat.v(end);
    [C_t_z, z, t, vr, Pr] = data_norm(dat,vmax);
catch
    vmax = input('vmax? \n');
    [C_t_z, z, t, vr, Pr] = data_norm(dat,vmax);
    realP = 0;
end

zz1 = [0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,0.95,0.98]; 

ci = @(ci0,vs,z,t)(ci0*(1-heaviside(z+vs*t)));
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
order = 3; frame = 15;
Ci = sgolayfilt(Ci,order,frame);


imax = find(1-abs(Ci/Ci(1)) > 0.001,1) - 1;
imin = find(abs(Ci/Ci(1)) > 0.001, 1, 'last' );

if isnan(imin)
    imin = length(t);
end

if imax < 5 
    imax = 10;
end

tvmin = t(imin);
Pvmin = abs(Ci(imin)/c0);
tvmax = t(imax);
Pvmax = 1-abs(Ci(imax)/c0);

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


[vs, P] = PVD_solve(solve,No);
length(P)

P = P/sum(P);
Nv = length(P);
stp = fix(length(t)/7);
Cmat = produce_data(P',vs,t(1:stp:end),z)';
Ctest = C_t_z(:,1:stp:end);
err(lmnop) = 1/length(Cmat(:))*sum((Cmat(:)-Ctest(:)).^2);

Cmatinerp = movmean(Cmat,20);
if pltC
figure
hold on
for col = 1:length(Cmat(1,:))
R = rand; G = rand; B = rand;
plot(Cmat(:,col),z,'-+','color',[R,G,B]); 
plot(C_t_z(:,(col-1)*stp+1),z,'color',[R,G,B]);
end
title('Concentration profiles')
xlim([0 c0*1.1]); ylim([min(z),max(z)])
xlabel('C [g/L]'); ylabel('z')

% legend('reconstructed data','input data');
set(gca,'FontSize',18)

if pltP
    figure
    plot(vs,P); hold on;
    plot(vr,Pr);
end
end

end

figure
plot(zz1,err);
toc