clear all;
close all;


sig1 = 1e-4 *[8 9 10 14];
var1 =[0.5 1.0 1.0 0.9]*1e-4;
Np1 = 1e8*[0.2 0.8 0.8 0.8];

Nv = 10000;


rhof = 1000;
rhop = 2500;
g = 9.81;
nu = 1e-6;

c0tot = 0.5;
z0 = 0.4;
z1 = 0.2;

a = z1/z0;


z = linspace(0,-0.5,400);

Nt = 200;

vv = linspace(40*(10^-4),1*10^(-5),200);
%t = z0./vv;%logspace(0,3,300);
t = 0:0.1:1000;
%t = linspace(1e-5,300,Nt);
%t = [linspace(1e-5,1e-4,50),linspace(1e-4,1e-3,50),linspace(1e-3,1e-2,50),linspace(1e-2,1e-1,50),linspace(1e-1,1,50),linspace(1,10,50),linspace(10,100,50)];

y = [];
for i = 1:length(sig1)
    y1 = randraw('normal', [sig1(i),var1(i)], Np1(i)); %Creates a random gaussian distribution of particles centered on diameter dp50 with variance var_dp
    y = [y ;y1];
end
[hist_part,vs] = hist(y,Nv); %Calculates the histogram of particle size with Nd bins, and returns the histogram and the mean particle size for each bin

figure
plot(vs,hist_part);
%plot(vs,(hist_part/max(hist_part))*max(P));

Np =sum(Np1);
ci0 = c0tot *hist_part./Np;

ci = @(ci0,vs,z,t)(ci0*(1-heaviside(z+vs*t)));


[~,iz0] = min((abs(z+z0)));
[~,iz1] = min((abs(z+z1)));


col = lines(20);

for k = 1:length(t)
    k
    C = 0;
    for i = 1:Nv
        C = ci(ci0(i),vs(i),z,t(k)) + C;
    end
    
    %figure(1)
    %plot(C,z); hold on;
    %C = C;%smooth(C,10);
    %plot(Cs,z);
    Cmat(k,:) = C;
   
end

     %   nme = ['Cdata_' num2str(k,'%05d') '.mat']

     Cdat.C =  Cmat;
     Cdat.t = t;
     Cdat.z = z;
     save(['created_data/Cdat1.mat'],'Cdat');
    %save(['created_data/' nme],'C');
    %save(['created_data/time.mat'],'t');
    %save(['created_data/depth.mat'],'z');

