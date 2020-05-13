clear 
close all
tic
addpath('./Toolboxes')

data2comp = 1:10;

for p = data2comp 
inA = ['./Marco/created_data/ML_Noisy/data_set_' num2str(p) '.mat'];
load(inA);

inM = ['../Model/VGG16/PVDs/ML_Noisy/PVDML_ds_'  num2str(p) '.txt'];
Pm = load(inM);

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


[vs, P] = PVD_solve(solve);

Pae{lmnop} = P/sum(P);
Nv = length(P);
stp = fix(length(t)/20);
Cmat = produce_data(P',vs,t(1:stp:end),z)';
Ctest = C_t_z(:,1:stp:end);
EAz(lmnop) = 1/length(Cmat(:))*sum((Cmat(:)-Ctest(:)).^2);
vsz{lmnop} = vs;
end

[val, id] = min(EAz);
Pa = Pae{id};
vs = vsz{id};
%%

v = linspace(0,1,100);

Pai = interp1(vs,Pa,v)./max(Pa);
Pai(isnan(Pai))=0;

Pri = interp1(vr,Pr,v)./max(Pr);
Pri(isnan(Pri))=0;

Nva = length(Pa);
Nvm = length(Pm);

ci0a = c0 * Pa;
ci0m = c0 * Pm/sum(Pm);

ci = @(ci0,vs,z,t)(ci0*(1-heaviside(z+vs*t)));

stp = fix(length(t)/40);
r = rem(length(t),stp);

erra = 0;
errm = 0;

for k = 1:stp:length(t)-r
    
    Ca = 0;
    for i = 1:Nva
        Ca = ci(ci0a(i),vs(i),z,t(k)) + Ca;
    end
    
     Cm = 0;
    for j = 1:Nvm
        Cm = ci(ci0m(j),v(j),z,t(k)) + Cm;
    end
    
    erra = sum(1./length(Ca)*(Ca'-C_t_z(:,k)).^2) + erra;
    errm = sum(1./length(Cm)*(Cm'-C_t_z(:,k)).^2) + errm;
end

ErraC(p-data2comp(1)+1) = erra;
ErrmC(p-data2comp(1)+1) = errm;

ErraP(p-data2comp(1)+1) = sum(1./100*(Pai - Pri).^2);
ErrmP(p-data2comp(1)+1) = sum(1./100*(Pm' - Pri).^2);

figure
plot(Pri,'--')
hold on
plot(Pai,'Linewidth',2)
plot(Pm,'Linewidth',2)
legend('Original PVD','Anlytic', 'ML')
set(gca,'FontSize',18)
end
%%
group = [1 * ones(size(ErraC));
         2 * ones(size(ErrmC))];
measure = [ErraC; ErrmC];
figure
boxplot(measure(:),group(:))
set(gca,'XTickLabel',{'Analytical','ML'})

group = [1 * ones(size(ErraP));
         2 * ones(size(ErrmP))];
measure = [ErraP; ErrmP];
figure
boxplot(measure(:),group(:))
set(gca,'XTickLabel',{'Analytical','ML'})

toc