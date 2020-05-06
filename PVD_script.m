clear; close all; tic
addpath('./Toolboxes')

pltC = 0; %display the concentration profiles
pltP = 0; % display the PVDs
pltE = 1; % display the 

in = 'infile';
data = importdata(in);


try
    t = data.t;
    z = data.z; % attention untis of z
    C_t_z = data.C;
    [m, ~] = size(C_t_z);
    if m ~= length(z)
        C_t_z = C_t_z';
    end
    if isfield(data,'P')
        Pr = data.P;
        vr = data.v;
        realP = 1;
    else
        realP = 0;
    end
catch
    
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
imax = find(1-abs(C_t_z(z1id,:)/c0) > 0.001, 1 )-1; 
imin = find(abs(C_t_z(z0id,:)/c0) > 0.001, 1, 'last' );

if isnan(imin)
    imin = length(t);
end

for i = 1:length(t)
    Ci(i)     = trapz(z(z1id:z0id),-C_t_z(z1id:z0id,i));
end

tvmin = t(imin);
Pvmin = abs(Ci(imin)/c0);
tvmax = t(imax);
Pvmax = 1-abs(Ci(imax)/c0);

K = ;
dv = K*z1*(1/t(imax+1)-1/t(imax));

solve.t =t;
solve.z = z;
solve.Ci = Ci;
solve.c0 = c0;
solve.z0 = z0;
solve.z1 = z1;
solve.tvmin = tvmin;
solve.tvmax = tvmax;
solve.dv = dv;
solve.Pvmin = Pvmin;
solve.Pvmax = Pvmax;


[vs, P] = PVD_solve(solve);
length(P)

P = P/sum(P);
Nv = length(P);
stp = fix(length(t)/20);
Cmat = produce_data(P',vs,t(1:stp:end),z)';
Ctest = C_t_z(:,1:stp:end);
err(lmnop) = 1/length(Cmat(:))*sum((Cmat(:)-Ctest(:)).^2);

if pltC
figure
plot(Cmat,z,'+'); hold on;
title('Concentration profiles')
xlim([0 c0*1.1]); ylim([min(z),max(z)])
plot(C_t_z(:,1:stp:end),z);
xlabel('C [g/L]'); ylabel('z [cm]')
legend('reconstructed data','input data');
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