clear 
close all
tic
%data_fabrication_DL

Disp = 1;

in = './TrainingData/Cdata_v10.mat';
data = load(in);

try 
    t = data.dat.t;
    z = data.dat.z; % attention untis of z
    C_t_z = data.dat.C;
    [m, ~] = size(C_t_z);
    if m ~= length(z)
        C_t_z = C_t_z';
    end
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

z1 = 0.95; 
z0 = 1;

z1dim = z1 * max(abs(z));
z0dim = z0 * max(abs(z));
sprintf('proposed z1 = %3.3f m\nproposed z0 = %3.3f m',z1dim,z0dim)

c0 = mean(C_t_z(:,1));
cN = mean(C_t_z(:,end));

C_v = C_t_z./c0.*(1-C_t_z./c0);

[~, z1id] = min( abs( z + z1dim) );
[~, z0id] = min( abs( z + z0dim) );

Ci = zeros(1,length(t)) ; %Ciwvl = Ci;
for i = 1:length(t)
    %wvlC=wden(C_t_z(:,i),'modwtsqtwolog','s','mln',4,'sym4');
     
    Ci(i)     = trapz(z(z1id:z0id),-C_t_z(z1id:z0id,i));
    %Ciwvl(i)  = trapz(z(z1id:z0id),-wvlC(z1id:z0id));
end

Ci = Ci./max(Ci); %Ciwvl = Ciwvl./max(Ciwvl);

dat  = struct('t',t,'z',z,'Ci',Ci,'c0',c0,'z0',z0dim,'z1'... 
                ,z1dim,'Disp',Disp,'C_t_z',C_t_z,'C_v',C_v); % Chose Ci or Ciwvl

[vs, P] = PVD_solve(dat);

if realP
    subplot(3,2,[4 6])
    plot(vr,Pr./max(Pr),'LineWidth',2,'LineStyle','--');
    legend('scatter','Fit','P original')
end

toc