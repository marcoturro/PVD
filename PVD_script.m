clear
%close all

Disp = 1
per = 20
tic

in='./TrainingData/Cdata_v10.mat';
data = load(in);

try 
    t = data.dat.t;
    z = data.dat.z; % attention exp data does not nned the 100 multiplier
    C_t_z = data.dat.C;
catch
end

ratio = per; % percentage
z1 = -floor(z(ceil(length(z)*ratio/100)));
z0 = -ceil(z(end - ceil(length(z)*ratio/100)));
c0 = mean(C_t_z(:,1));
cN = mean(C_t_z(:,end));
P0 = cN/c0;

sprintf('proposed z0 = %d cm\nproposed z1 = %d cm',z0,z1)
%pause()
[~, z1id] = min( abs( z + z1) );
[~, z0id] = min( abs( z + z0) );


for i = 1:length(t)
    wvlC=wden(C_t_z(:,i),'modwtsqtwolog','s','mln',4,'sym4');
     
    Ci     = [Ci     sum(-C_t_z(z1id:z0id,i)'.*dz(z1id:z0id))];
    Ciwvl  = [Ciwvl  sum(-wvlC(z1id:z0id).*dz(z1id:z0id))];
end

Vct  = {t,z,Ci,P0,c0,z0,z1,z1id,z0id,C_t_z}; % Chose Ci or Ciwvl

[vs, P, ERR] = PVD_solve(Vct,Disp);

figure
plot(vs,P)

toc