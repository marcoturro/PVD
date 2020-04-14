%close all
%clear

%make into a loop for higher number of datasets
nbOfSets = 1;
for s = 1:nbOfSets
    
%%
[vad, Pv] = randP(randi([90 140],1)); %to test
vs = 0.04*randi([30 50],1)*vad; % for realism


c0tot = 0.5;
t = linspace(0,randi([1000 1400],1),randi([8000 14000],1)); %to be twitched
ci0 = c0tot * Pv;
ci = @(ci0,vs,z,t)(ci0*(1-heaviside(z+vs*t)));
C = [];
z = 0:-0.01*randi([5 10],1):-10*randi([3 5],1);
%%
C = zeros(length(z),length(t));
for k = 1:length(t)
    t(k)/t(end)*100
    Ci = 0;
    for i = 1:length(Pv)
        Ci = ci(ci0(i),vs(i),z,t(k)) + Ci;
    end
%     plot(Ci,z)
%     xlim([0 c0tot]);
%     ylim([min(z) max(z)]);
%     pause(0.1)
    C(:,k) = Ci';
    
end


fileID_Cdata= ['./TrainingData/Cdata_v' num2str(s) '.mat'];
dat = struct('C',C,'t',t,'z',z,'P',Pv,'v',vad);
save(fileID_Cdata,'dat')

end
