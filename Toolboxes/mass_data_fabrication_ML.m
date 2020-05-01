clear all;
close all;
addpath 'E:\Raphael\MIT\sediment_model\GIT\PVD\Raphael'

if ~exist('E:\Raphael\MIT\sediment_model\ML_TRAINING_DATA','dir');
mkdir 'E:\Raphael\MIT\sediment_model\ML_TRAINING_DATA';
else
    nnn = dir('E:\Raphael\MIT\sediment_model\ML_TRAINING_DATA\*.mat');
    nnn = numel(nnn);
end


Np = 200;

vmax = 1;
vmin = 10^(-4);

Nz = 100;

z0 = 1;

a2 = 4000;
b2 = 12000;

a3 = 0.1;
b3 = 0.4;

apol = 2;
bpol = 5;
%Nsamples = 1e7;



Ndat = 100;

z = linspace(0,-1,Nz);
cnt_max = 10000;
cnt = 0;
while cnt < cnt_max
    cnt = cnt+1
vmin = (b3-a3)*rand+a3;
tend = 1/vmin;
Nt =(b2-a2)*rand + a2; 
dt = tend/Nt;

t = 0:dt:tend;
v = linspace(vmin,1,Np);

Npol = randi([apol,bpol]);
xx = [vmin rand(1,Npol) 1];
yy = [0 rand(1,Npol) 0];
p = polyfit(xx,yy,Npol+1);


P = polyval(p,v).^2.*(1-0.2*rand(1,Np));
P = P/sum(P);

% fprintf('vmin is %3.3f\n',vmin);
% fprintf('dt is %0.7f\n',dt);
% fprintf('tmax is %3.3f\n',tend);


Cmat = produce_data(P,v,t,z);

dat.C = Cmat;
dat.z = z;
dat.t = t;
dat.P = P;
dat.v = v;

save(['E:\Raphael\MIT\sediment_model\ML_TRAINING_DATA\data_set_' num2str(nnn+cnt) '.mat'],'dat');


end


%ci = produce_data(P,t,z);