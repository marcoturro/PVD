close all
clear
addpath('./Toolboxes')
 

sets = 1:1;

tic
for s = sets
   clearvars -except s sets
%%
ref = randi([500 550],1);
[vad, Pv] = randP(ref); 
vs = linspace(0,0.0005*randi([20,30],1),length(vad)+1); % for realism
vs = vs(2:end);
vmax = randi([9 50])/1000;
vs = vad*vmax;

c0tot = randi([1 5],1)/10;
t = linspace(0,randi([8000 10000],1),randi([1000 1200],1)); %to be twitched
ci0 = c0tot * Pv;
z = 0:-0.01*randi([8 10],1):-10*randi([2 5],1);
z = z./100;

ci = @(ci0,vs,z,t)(ci0*(1-heaviside(z+vs*t)));
ci2 =  @(ci0,vs,z,t)(ci0*sigmoid(z,vs*t));
C = [];

%%
C = zeros(length(z),length(t));
for k = 1:length(t)
    %t(k)/t(end)*100
    Ci = 0;
    for i = 1:length(Pv)
        Ci = ci(ci0(i),vs(i),z,t(k)) + Ci;
    end
    Ci = movmean(Ci,7);
    f_1 = 5000;
    a_1 = 600;
    
    f_2 = 0.01;
    a_2 = 1/100;
    g_2 = 1.5*max(abs(z));
    v_2 = f_2*g_2;
    
    f_3 = 0.05;
    a_3 = 1/250;
    g_3 = 0.5*max(abs(z));
    v_3 = f_2*g_2;
    
   Cin = awgn(Ci,48) + Ci.*awgn(Ci,60)/5; % + Ci.*awgn(sin((z-v_2*t(k))*2*pi/g_2)*a_2,65) + Ci.*awgn(sin((z-v_3*t(k))*2*pi/g_3)*a_3,65);
   
%        plot(Ci,z)
%        xlim([0 1.3*c0tot]);
%        ylim([min(z) max(z)]);
%        pause(0.5)
%     
    Cn(:,k) = Cin';
    C(:,k) = Ci';
    
end

%%
fileID_Cdata= ['./Marco/created_data/data_set_' num2str(s) '.mat'];
dat = struct('C',C,'t',t,'z',z,'P',Pv,'v',vs);
save(fileID_Cdata,'dat')

% fileID_Cdata= ['../ML_TRAINING_DATA/ML_MT/Noisy_b2/data_set_' num2str(s) '.mat'];
% dat = struct('C',Cn,'t',t,'z',z,'P',Pv,'v',vs);
% save(fileID_Cdata,'dat')

sprintf([num2str(s) ' out of ' num2str(sets(end))])

end
toc