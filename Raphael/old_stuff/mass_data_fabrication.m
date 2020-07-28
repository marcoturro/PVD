clear all;
close all;


%Nt = 200;
Np = 200;

vmax = 1;
vmin = 10^(-4);

Nz = 400;

z0 = 1;

a2 = 4000;
b2 = 20000;

a3 = 0.1;
b3 = 0.4;

Nsamples = 1e7;
a = 0.5;
b = 4;
omega = 1;

Ndat = 100;

for i = 1:Ndat
    i
    %z1 =(b3-a3)*rand+a3;
    z = linspace(0,-1,Nz);

% gamma = 10^((b2-a2).*rand+a2);
%v_for_t = linspace(1,1/gamma,Nt);
%t = 1./v_for_t;
vmin = (b3-a3)*rand+a3;
tend = 1/vmin;
Nt =(b2-a2)*rand + a2; 
dt = tend/Nt;

t = 0:dt:tend;
m = (b-a).*rand + a;

[x] = randraw('nakagami', [m, omega], Nsamples);
[P,wrk] = hist(x,Np);

fprintf('vmin is %3.3f\n',vmin);
fprintf('dt is %0.7f\n',dt);
fprintf('tmax is %3.3f\n',tend);

P = P/Nsamples;
P = fliplr(P);
wrk = wrk /max(wrk);
figure(1)
plot(wrk,P); hold on;

Cmat = produce_data(P,t,z);
figure(2)
plot(t,sum(Cmat,2)); hold on;

dat.C = Cmat;
dat.z = z;
dat.t = t;
dat.P = P;
dat.v = wrk;

save(['created_data/data_set_' num2str(i) '.mat'],'dat');


end


%ci = produce_data(P,t,z);