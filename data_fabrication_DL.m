
%make into a loop for higher number of datasets
nbOfSets = 500;
tic
for s = 1:nbOfSets
   clearvars -except s nbOfSets 
%%
[vad, Pv] = randP(randi([90 140],1)); %to test
vs = linspace(0,0.000005*randi([20,30],1),length(vad)+1); % for realism
vs = vs(2:end);

close

c0tot = 0.5;
t = linspace(0,randi([8000 10000],1),randi([16000 20000],1)); %to be twitched
ci0 = c0tot * Pv;
z = 0:-0.01*randi([4 8],1):-10*randi([3 5],1);
z = z./100;

ci = @(ci0,vs,z,t)(ci0*(1-heaviside(z+vs*t)));
ci2 =  @(ci0,vs,z,t)(ci0*heaviside(z+vs*t));
C = [];

%%
C = zeros(length(z),length(t));
for k = 1:length(t)
    %t(k)/t(end)*100
    Ci = 0;
    for i = 1:length(Pv)
        Ci = ci(ci0(i),vs(i),z,t(k)) + Ci;
    end
%     plot(Ci,z)
%     xlim([0 2*c0tot]);
%     ylim([min(z) max(z)]);
%     pause(0.1)
    C(:,k) = Ci';
    
end


fileID_Cdata= ['./TrainingData/data_set_' num2str(s) '.mat'];
dat = struct('C',C,'t',t,'z',z,'P',Pv,'v',vs);
save(fileID_Cdata,'dat')

sprintf(num2str(s))

end
toc