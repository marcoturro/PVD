clear all;
close all;

% for n = 1:42
% fle = ['Raphael/created_data/data_set_' num2str(n) '.mat'];
% 
% dat = importdata(fle);
% 
% figure(1)
% plot(dat.v,dat.P); hold on;
% xlabel('$v$');
% ylabel('$P$');
% end


fold = 'E:\Raphael\MIT\sediment_model\ML_TRAINING_DATA';


for n = 1750:1767
dat = importdata([fold '\data_set_' num2str(n) '.mat']);

stp = 200;

figure(1)
plot(dat.v,dat.P); hold on;
end