clear
% close all

dict = {'Bottle1-1'};
dict = {'Marco05.1.7.20','Marco05.2.7.20','Marco05.3.7.20','Marco05.4.7.20','Marco05.6.7.20','Marco05.7.7.20'};
dict = {'05_sat_180','05_sat_240'}; 
dict = {'Marco05.3.7.20','Marco05.4.7.20'};
folder = './exp_data/';

figure; subplot(1,2,1); hold on
col = lines(length(dict))
all_marks = {'o','+','*','.','x','s','d','^','v','>','<','p','h'};

for i = 1:length(dict)
    load([folder dict{i} '.mat']);

    plot(dat.Imean_i(160:20:end,1:6),dat.z(160:20:end),'color',col(i,:)) %,'Marker',all_marks{i})
    Io(i) = dat.I0;
    I1(i) = mean(dat.Imean(:,3));
end
ylabel('z [m]')
xlabel('I [-]')
set(gca,'Fontsize',14)

%%
subplot(1,2,2); hold on
for i = 1:length(dict)
    load([folder dict{i} '.mat']);
    plot(dat.t(1:4:end)/3600,mean(dat.Imean(200:end,1:4:end),1),'color',col(i,:))
end
xlabel('time [h]')
ylabel('I [-]')
set(gca,'Fontsize',14)

