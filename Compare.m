clear 
close all
tic
addpath('./Toolboxes')

data2comp = 6001:8000;

for p = data2comp 
inA = ['../ML_TRAINING_DATA/ML_RO/data_set_' num2str(p) '.mat'];
load(inA);

inM = ['../Model/CNN_time_Model/PVDs/PVDML_ds_'  num2str(p) '.txt'];
Pm = load(inM);

try 
    t = dat.t;
    z = dat.z; % attention untis of z
    C_t_z = dat.C;
    [m, ~] = size(C_t_z);
    if m ~= length(z)
        C_t_z = C_t_z';
    end
    if isfield(dat,'P')
        Pr = dat.P;
        vr = dat.v;
        vr = dat.v/max(dat.v);
        realP = 1;
    else
        realP = 0;
    end
catch
    
end

z1 = 0.95; 
z0 = 1;

c0 = mean(C_t_z(:,1));
cN = mean(C_t_z(:,end));

C_v = C_t_z./c0.*(1-C_t_z./c0);

z1dim = z1 * max(abs(z));
z0dim = z0 * max(abs(z));

[~, z1id] = min( abs( z + z1dim) );
[~, z0id] = min( abs( z + z0dim) );

Ci = zeros(1,length(t)) ;

for i = 1:length(t)
     
    Ci(i)     = trapz(z(z1id:z0id),-C_t_z(z1id:z0id,i));

end

Ci = Ci./max(Ci);

clstr  = struct('t',t,'z',z,'Ci',Ci,'c0',c0,'z0',z0dim,'z1'... 
                ,z1dim,'Disp',0,'C_t_z',C_t_z,'C_v',C_v);

[vs, Pa, ERR] = PVD_solve(clstr);

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

end

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