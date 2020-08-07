clear
close all

dict = {'Bottle1-1'};
dict = {'Marco05.1.7.20','Marco05.2.7.20','Marco05.3.7.20','Marco05.4.7.20','Marco05.6.7.20','Marco05.7.7.20'};
dict = {'05_sat_180'}; 
dict = {'05_rotated','05_sat_180'};  
folder = 'D:\Sediment_exp\Measurements\';

figure; subplot(1,2,1); hold on
col = lines(length(dict))
all_marks = {'o','+','*','.','x','s','d','^','v','>','<','p','h'};

for i = 1:length(dict)
    load([folder dict{i} '\' dict{i} '.mat']);

    plot(dat.Imean_i(1:20:end,1:24),dat.z(1:20:end),'color',col(i,:)) %,'Marker',all_marks{i})
    Io(i) = dat.I0;
    I1(i) = mean(dat.Imean(:,3));
end
ylabel('z [m]')
xlabel('I [-]')
set(gca,'Fontsize',14)


subplot(1,2,2); hold on
for i = 1:length(dict)
    load([folder dict{i} '\' dict{i} '.mat']);
    plot(dat.t(1:4:end)/3600,mean(dat.Imean(200:end,1:4:end),1),'color',col(i,:))
end
xlabel('time [h]')
ylabel('I [-]')
set(gca,'Fontsize',14)

name = dict{1}
suptitle(['C = ' name(6:7) ' g/L'])