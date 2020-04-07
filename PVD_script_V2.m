clear
%close all

D=4
per = 30
tic

%in = ['./Cdata/Cdata_v' num2str(D) '.mat']; %input('data file name? \nEx. Cdata.mat \n','s');
 in='./Cdata/Cdata_v1.mat'
data = load(in);

try 
    t = data.dat.t;
    z = data.dat.z*100;
    C_t_z = data.dat.C;
catch
end



figure; subplot(2,2,1); plot(t); title('data time [s]'); subplot(2,2,3); dt = mean(diff(t)); plot(diff(t)); subplot(2,2,2); plot(z); title('space [cm]'); subplot(2,2,4); dz = diff(z); plot(diff(z))

ratio = per; % percentage
z1 = -floor(z(ceil(length(z)*ratio/100)));
z0 = -ceil(z(end - ceil(length(z)*ratio/100)));
c0 = mean(C_t_z(:,1));
cN = mean(C_t_z(:,end));
P0 = cN/c0;

sprintf('proposed z0 = %d cm\nproposed z1 = %d cm',z0,z1)
%pause()
[~, z1id] = min( abs( z + z1) );
[~, z0id] = min( abs( z + z0) );

Ciwvl= []; Cimavg = []; Ci = [];

for i = 1:length(t)
    
    wvlC=wden(C_t_z(:,i),'modwtsqtwolog','s','mln',4,'sym4');
    mavgC=movmean(C_t_z(:,i),length(z)*0.05)';
    
    t(i);
    clc
    
    Ci     = [Ci     sum(-C_t_z(z1id:z0id,i)'.*dz(z1id:z0id))];
    Ciwvl  = [Ciwvl  sum(-wvlC(z1id:z0id).*dz(z1id:z0id))];
    Cimavg = [Cimavg sum(-mavgC(z1id:z0id).*dz(z1id:z0id))];
end

figure
hold on
plot(t,Ci)
plot(t,Ciwvl)
plot(t,Cimavg); hold off
legend('Ci','Ci Wavelet smoothing','Ci moving average')
P0=0;

Vct  = {t,z,Cimavg,P0,c0,z0,z1,C_t_z};
%%
rng default % for reproducibility]
objconstr = @(x)PVD_Opt(x,Vct);
IntCon = 1:2; lb = [1 1]; ub = [30 600];
opts = optimoptions(@ga, ...
                             'PopulationSize', 15, ...
                             'MaxGenerations', 10, ...
                             'EliteCount', 2, ...
                             'FunctionTolerance', 1e-3); %...
                             %'PlotFcn', @gaplotbestf);
[xopt,fval,exitflag,Output] = ga(objconstr,2,[],[],[],[],...
lb,ub,[],IntCon,opts);

xopt
fval
toc

%% Visual run
% close all; PVD_Opt(xopt,Vct,1)

figure
plot(t,Ci)
title([in(1:end-7) ' ' in(end-5:end-4)])
yl = ylim;
ylim(yl)
hold on
plot(t*0+xopt(2),t)
legend('Ci')
[~, index] = min(abs(t-xopt(2)));
mLoss = (1-Ci(index)/Ci(1))*100;
lossRate = mLoss/xopt(2);
set(gca,'FontSize',18)
txt{1} = cellstr(['  tSwitch linear-log = ' num2str(xopt(2)) ' s']);
txt{2} = cellstr(['  Mass Loss = ' num2str(mLoss) ' %']);
txt{3} = cellstr(['  Linear refinement = ' num2str(xopt(1)) ' steps']);
txt{4} = cellstr(['  Mass Loss Rate = ' num2str(lossRate) ' %.s-1']);
txt{5} = cellstr(['  MS Error = ' num2str(fval) ' [ - ]']);
txt{6} = cellstr(['  z1 = ' num2str(z1) ' , z0 = ' num2str(z0) ' [cm]']);
text(xopt(2), yl(2)*0.95,txt{1}, 'Fontsize', 14);
text(xopt(2), yl(2)*0.9,txt{2}, 'Fontsize', 14);
text(xopt(2), yl(2)*0.85,txt{3}, 'Fontsize', 14);
text(xopt(2), yl(2)*0.8,txt{4}, 'Fontsize', 14);
text(xopt(2), yl(2)*0.75,txt{5}, 'Fontsize', 14);
text(xopt(2), yl(2)*0.7,txt{6}, 'Fontsize', 14);
