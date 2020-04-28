clear
close all

m = length(dat.z)+1 %normally length(dat.z) = 100
n = 100

for j = 1:300
load(['../ML_TRAINING_DATA/data_set_' num2str(j) '.mat']);
if length(dat.C(:,1)) ~= length(dat.z)
    dat.C = dat.C';
end

C = zeros(m,n);
C(1,:) = interp1(linspace(0,1,length(dat.t)),dat.t,linspace(0,1,n))/max(dat.t);
for i = 2:length(dat.z+1)
    sel = dat.C(i,dat.C(i,:)>0.0001 & dat.C(i,:)<0.9999);
    C(i,:) = interp1(linspace(0,1,length(sel)),sel,linspace(0,1,n))';
end
dat.C = C;
surf(dat.C)
save(['./ML_Preprocessed_CNN_100_time/data_set_' num2str(j) '.mat'],'dat')
end


