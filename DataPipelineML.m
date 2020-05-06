clear
close all
kernel = 3
Pad = 1
Isize = 100
m = Isize + kernel + 2*Pad
n = Isize + 2*Pad
type = 'linear'
Data_Pool = 'ML_TRAINING_DATA/ML_RO'


for DatN = 6000:8000
    
load(['../' Data_Pool '/data_set_' num2str(DatN) '.mat'])

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
            C(i,1+Pad:end-Pad) = interp1(linspace(0,dat.t(ti),ti),dat.C(i-kernel - Pad,1:ti),dx)';
        end
      
end

dat = rmfield(dat,'C');
dat.C = C;
save(['../ML_TRAINING_DATA/ML_Preprocessed_CNN_time_norm/data_set_' num2str(DatN) '.mat'],'dat')
end


        
%%%%% test for varying long dt %%%%
%         for k = 1:m-1
%             ti(k) = length(dat.C(k,dat.C(k,:)>0.002)); %vector
%             dx(k,:) = 1./linspace(1/0.1,1/dat.t(ti(k)),n); %nonlinear in time
%         end
%         C = zeros(m,n);
%         C(1,:) = interp1(linspace(0,1,length(dat.t(1:ti(end)))),dat.t(1:ti(end)),linspace(0,1,n))/10;
% 
%         for k = 2:m
%             try
%             C(k,:) = interp1(linspace(0,dat.t(end),ti(k-1)),dat.C(k-1,1:ti(k-1)),dx(k-1,:))';
%             catch
%             C(k,1) = dat.C(k,1);
%             end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%     case 'nonlinear'
%        ti = length(dat.C(end,dat.C(end,:)>0.002)); %scalar
%         C = zeros(m,n);
%         C(1,:) = interp1(linspace(0,1,length(dat.t(1:ti))),dat.t(1:ti),linspace(0,1,n))/10;
%         dx = [0 1./linspace(1/0.1,1/dat.t(ti),n-1)]; %nonlinear 
%         for i = 2:m
%             C(i,:) = interp1(linspace(0,dat.t(end),ti),dat.C(i-1,1:ti),dx)';
%         end
