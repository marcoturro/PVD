clear all;
close all;
%close all
tic
%data_fabrication_DL
col = lines(100);
ci = @(ci0,vs,z,t)(ci0*(1-heaviside(z+vs*t)));

for itst = 1:42
in=['./Raphael/created_data/data_set_' num2str(itst) '.mat'];%./TrainingData/Cdata_v10.mat';
%in = ['./TrainingData/Cdata_v10.mat'];
data = load(in);

try 
    t = data.dat.t;
    z = data.dat.z; % attention untis of z
    C_t_z = data.dat.C;
    [m, ~] = size(C_t_z);
    if m ~= length(z)
        C_t_z = C_t_z';
    end
    %maxz = max(abs(z));
    %z = z./maxz;
    %maxt = max(t);
    %t = t./max(t);
    if isfield(data.dat,'P')
        Pr = data.dat.P;
        vr = data.dat.v;
        %vr = vr./(maxz/maxt);
        realP = 1;
    else
        realP = 0;
    end
catch
    
end

z1 = 0.9; %-(z(floor(length(z)*(ratio-offset)/100)));
z0 = 0.95;%-(z(end - floor(length(z)*(ratio+offset)/100)));

z1dim = z1 * max(abs(z));
z0dim = z0 * max(abs(z));
sprintf('proposed z1 = %3.3f m\nproposed z0 = %3.3f m',z1dim,z0dim)

c0 = mean(C_t_z(:,1));
cN = mean(C_t_z(:,end));

[~, z1id] = min( abs( z + z1dim) );
[~, z0id] = min( abs( z + z0dim) );

Ci = zeros(1,length(t)) ; Ciwvl = Ci;
for i = 1:length(t)
    %wvlC=wden(C_t_z(:,i),'modwtsqtwolog','s','mln',4,'sym4');
     
    Ci(i)     = trapz(z(z1id:z0id),-C_t_z(z1id:z0id,i));
    %Ciwvl(i)  = sum(-wvlC(z1id:z0id).*dz(z1id:z0id));
end

Ci = Ci./max(Ci); %Ciwvl = Ciwvl./max(Ciwvl);

dat  = struct('t',t,'z',z,'Ci',Ci,'c0',c0,'z0',z0dim,'z1'... 
                ,z1dim,'z1id',z1id,'z0id',z0id,'C_t_z',C_t_z); % Chose Ci or Ciwvl

[vs, P] = PVD_solve(dat);

if realP
    %figure
    %plot(vr,Pr./max(Pr),'LineWidth',2,'LineStyle','--'); hold on;
    %plot(vs,P/max(P));
    %legend('scatter','P original')
end

 t_test = linspace(t(1),t(end),50);
 figure
 for k = 1:length(t_test)
     
     C_original = 0;
     for i = 1:length(vr)
         C_original = ci(Pr(i)/sum(Pr),vr(i),z,t_test(k)) + C_original;
     end
     
          C_reconstructed = 0;
     for i = 1:length(vs)
         C_reconstructed = ci(P(i)/sum(P),vs(i),z,t_test(k)) + C_reconstructed;
     end
     err(k) = 1/length(z)*sum((C_original-C_reconstructed).^2);
     plot(C_original,z,'color',col(k,:)); hold on;
     plot(C_reconstructed,z,'+','color',col(k,:));

 end
ERR(itst) = max(err);

end

figure
plot(ERR,'o');
ylabel('$\epsilon$');
xlabel('id');

toc