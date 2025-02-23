clear
% close all

dict = {'Bottle1-1'};
dict = {'CCFZ_05.7.7.20','CCFZ_05_sat_900'};
folder = './exp_data/';

figure; hold on
col = lines(length(dict))
all_marks = {'o','+','*','.','x','s','d','^','v','>','<','p','h'};

for i = 1:length(dict)
    subplot(1,2,i)
    load([folder dict{i} '.mat']);
    if i == 1
        offs = 0;
    else
        offs = .008;
    end

    plot(dat.Imean_i(160:20:end,1:12),dat.z(160:20:end)+offs,'color',col(i,:)) %,'Marker',all_marks{i})
    Io(i) = dat.I0;
    I1(i) = mean(dat.Imean(:,3));
end
ylabel('z [m]')
xlabel('I [-]')
set(gca,'Fontsize',14)

%%
figure
subplot(1,2,2); hold on
for i = 1:length(dict)
    load([folder dict{i} '.mat']);
    plot(dat.t(1:4:end)/3600,mean(dat.Imean(200:end,1:4:end),1),'color',col(i,:))
end
xlabel('time [h]')
ylabel('I [-]')
set(gca,'Fontsize',14)

