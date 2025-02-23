clear
close all
addpath('./Toolboxes')

kernel = 3
Pad = 1
Isize = 100
m = Isize + kernel + 2*Pad
n = Isize + 2*Pad
type = 'linear'
Data_Pool = './exp_data/'
Save = './SiCML/';
exp = {'SiCF500_05.1.8.14','SiCF500_05.2.8.14','SiCF800_05.1.8.14','SiCF800_05.2.8.14','SiCF1000_05.1.8.14','SiCF1000_05.2.8.14','SiC33_05.1.8.16','SiC33_05.2.8.16','SiC50_05.1.8.17','SiC50_05.2.8.17'}
for DatN = 1:length(exp)
DatN

%load([Data_Pool '/data_set_' num2str(DatN) '.mat'])
load([Data_Pool exp{DatN} '.mat'])

try
    vmax = dat.v(end);
    [dat.C, dat.z, dat.t, dat.v, dat.P] = data_norm(dat,vmax);
catch
    vmax = input('vmax? \n');
    [dat.C, dat.z, dat.t, dat.v, dat.P] = data_norm(dat,vmax);
end

switch type
    case 'linear'
        ti = length(dat.C(end,dat.C(end,:)>0.01)); %scalar
        C = zeros(m,n);
        C(1+Pad:kernel+Pad,1+Pad:end-Pad) = ones(kernel,1)*interp1(linspace(0,1,length(dat.t(1:ti))),dat.t(1:ti),linspace(0,1,Isize))/10;
        dx = linspace(0,dat.t(ti),Isize); %linear 
        for i = kernel + Pad +1 : kernel + Pad + Isize
            C(i,1+Pad:end-Pad) = interp1(dat.t(1:ti),dat.C(i - kernel - Pad,1:ti),dx)';
        end
      
end

dat = rmfield(dat,'C');
dat.C = C;
save([Save '/data_set_' num2str(DatN) '.mat'],'dat')
end
