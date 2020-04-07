clear
%close all

Disp = 1
ratio = 40
offset = -35 % (positive is upwards ie. opposite of gravity accel)
tic

in='./TrainingData/Cdata_v10.mat';
data = load(in);

try 
    t = data.dat.t;
    z = data.dat.z; % attention untis of z
    C_t_z = data.dat.C;
%     maxz = max(abs(z));
%     z = z./maxz;
%     maxt = max(t);
%     t = t./max(t);
catch
end

z1 = -floor(z(ceil(length(z)*(ratio-offset)/100)));
z0 = -ceil(z(end - ceil(length(z)*(ratio+offset)/100)));
c0 = mean(C_t_z(:,1));
cN = mean(C_t_z(:,end));
P0 = cN/c0;
dz = diff(z);

sprintf('proposed z0 = %d cm\nproposed z1 = %d cm',z0,z1)

[~, z1id] = min( abs( z + z1) );
[~, z0id] = min( abs( z + z0) );

Ci = []; Ciwvl = [];
for i = 1:length(t)
    wvlC=wden(C_t_z(:,i),'modwtsqtwolog','s','mln',4,'sym4');
     
    Ci     = [Ci     sum(-C_t_z(z1id:z0id,i)'.*dz(z1id:z0id))];
    Ciwvl  = [Ciwvl  sum(-wvlC(z1id:z0id).*dz(z1id:z0id))];
end

dat  = struct('t',t,'z',z,'Ci',Ci,'P0',P0,'c0',c0,'z0',z0,'z1'... 
                ,z1,'z1id',z1id,'z0id',z0id,'C_t_z',C_t_z); % Chose Ci or Ciwvl

[vs, P, ERR] = PVD_solve(dat,Disp);

toc