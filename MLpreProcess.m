clear
close all
m = 101
n = 100


for j = 320:6000
load(['../ML_TRAINING_DATA/data_set_' num2str(j) '.mat']);
if length(dat.C(:,1)) ~= length(dat.z)
    dat.C = dat.C';
end

if length(dat.z) ~= m - 1
    error('the data is not normalised in space (100 points [0 ; 1])')
end

ti = length(dat.C(end,dat.C(end,:)>0.01));
C = zeros(m,n);
C(1,:) = interp1(linspace(0,1,length(dat.t(1:ti))),dat.t(1:ti),linspace(0,1,n))/10;

for i = 2:m
    C(i,:) = interp1(linspace(0,1,ti),dat.C(i-1,1:ti),linspace(0,1,n))';
end
dat = rmfield(dat,'C');
dat.C = C;
save(['./ML_TRAINING_DATA/ML_Preprocessed_CNN_' num2str(m) '/data_set_' num2str(j) '.mat'],'dat')
end