function [v, P, ERR] = PVD_solve(dat,Disp)
%%
clearvars -except dat Disp

f_t = dat.t;
f_z = dat.z;
f_Ci = dat.Ci;
f_P0 = dat.P0;
f_c0 = dat.c0;
f_z0 = dat.z0;
f_z1 = dat.z1;
f_z1id = dat.z1id;
f_z0id = dat.z0id;
f_C_t_z = dat.C_t_z;

 if ~exist('Disp','var')
     % third parameter does not exist, so default it to something
      Disp = 0;
 end

order = 3; frame = 13;
f_Ci(2:end) = sgolayfilt(f_Ci(2:end),order,frame);
f_Ci(2:end)= wdenoise(f_Ci(2:end),4);

for i = 1:length(f_Ci)
    if  abs((f_Ci(i)-f_Ci(1)))/f_Ci(1)*100 >= 0.01
        tvmax = f_t(i);
        dv = 15*f_z1*(1/f_t(i+1)-1/f_t(i));
        break
    end
    
end

for j = length(f_Ci):-1:i
    if  abs((f_Ci(end)-f_Ci(j))) >= 0.01
        tvmin = f_t(j);
        break
    end
end
 
vmax= f_z1/tvmax ;
vmin = f_z1/tvmin;
vt = vmax:dv:vmin;
ti = f_z1./vt;

tint = zeros(1,length(ti));
Id = tint;

for i = 1:length(ti)  
    [~, idx] = min(abs(f_t-ti(i)));
    Id(i) = idx;
    tint(i) = f_t(idx);
end

dz = diff(f_z);
Ci = zeros(1,length(tint));
for i = 1:length(tint)
    Ci(i) = sum(-f_C_t_z(f_z1id:f_z0id,Id(i))'.*dz(f_z1id:f_z0id));
end

%Ci_i = interp1(f_t,f_Ci,ti);
Ci_i = Ci;
ti = tint;

try
order = 3; frame = 13;
Ci_i = sgolayfilt(Ci_i,order,frame);
catch
order = 3; frame = 2*floor(length(ti))/2-1;
Ci_i = sgolayfilt(Ci_i,order,frame);
end

dddtCint = my_2FD_non_uniform(ti,Ci_i);
dddtCint(dddtCint<0)=0;
[v, P] = P_RHS_solve(f_z0,f_z1,ti,dddtCint,f_P0);

Nv = length(P);
ci0 = f_c0 * P;
ci = @(ci0,vs,z,t)(ci0*(1-heaviside(z+vs*t)));
ERR = 0;

if Disp == 1
figure('units','normalized','outerposition',[0 0 1 1])
end

div = 40;
r = rem(length(f_t),div);
stp = (length(f_t)-r)/div;

for k = 1:stp:(length(f_t)-r)
    
    C = 0;
    for i = 1:Nv
        C = ci(ci0(i),v(i),f_z,f_t(k)) + C;
    end
       err = sum((C-movmean(f_C_t_z(:,k),length(f_z)*0.05)').^2); % OVERALL ERROR
      
       ERR = ERR + err;
       
       if Disp == 1
       [~, boundD] = min(abs(f_C_t_z(:,k)-f_c0*0.2)) ; [~, boundU] = min(abs(f_C_t_z(:,k)-f_c0*0.8));
       smallErr = sum((C(1:boundD)-movmean(f_C_t_z(1:boundD,k),length(f_z(1:boundD))*0.05)').^2);
       bigErr = sum((C(boundU:end)-movmean(f_C_t_z(boundU:end,k),length(f_z(boundU:end))*0.05)').^2);
       sb1 = subplot(3,2,[1 3 5]) ;
       plot(movmean(C,3),f_z,'--')
       hold on
       title('Concentration profiles')
       xlim([0 f_c0*1.1]); ylim([min(f_z),max(f_z)])
       plot(movmean(f_C_t_z(:,k),length(f_z)*0.05)',f_z)
       scatter([f_C_t_z(boundD,k) f_C_t_z(boundU,k)],[f_z(boundD) f_z(boundU)])
       legend('reconstructed data','input data','Measure points')
       set(gca,'FontSize',18)
       sb2 = subplot(3,2,2);
       scatter(f_t(k),smallErr)
       xlim([0 f_t(end)])
       title('mean square difference for small particles')
       legend(['time remaining : ' sprintf('%.0f',f_t(end)-f_t(k)) ' [s]'])
       set(gca,'FontSize',18)
       sb3 = subplot(3,2,4);
       scatter(f_t(k),bigErr)
       xlim([0 f_t(end)])
       title('mean square difference for large particles')
       %legend(['time remaining : ' sprintf('%.0f',f_t(end)-f_t(k)) ' [s]'])
       set(gca,'FontSize',18)
       sb4 = subplot(3,2,6);
       plot(v,P,'LineWidth', 3)
       title('Particle Velocity Distribution')
       set(gca,'FontSize',18)
%        scatter(f_t(k),err)
%        xlim([0 f_t(end)])
%        title('mean square difference for all particles')
%        legend(['time remaining : ' sprintf('%.0f',f_t(end)-f_t(k)) ' [s]'])
       hold(sb2,'on') ; hold(sb3,'on') ; hold(sb4,'on')
       pause (0.1)
       end
       if err <= 0.0001
        break
       end
            
end




end

