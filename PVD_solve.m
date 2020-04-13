function [v, P] = PVD_solve(dat)
%%
clearvars -except dat Disp

t = dat.t;
z = dat.z;
Ci = dat.Ci;
c0 = dat.c0;
z0 = dat.z0;
z1 = dat.z1;
C_t_z = dat.C_t_z;
C_v = dat.C_v;
Disp = dat.Disp;

 if ~exist('Disp','var')
     % third parameter does not exist, so default it to something
      Disp = 0;
 end

% order = 3; frame = 13;
% f_Ci(2:end) = sgolayfilt(f_Ci(2:end),order,frame);
% f_Ci(2:end)= wdenoise(f_Ci(2:end),4);

for i = 1:length(Ci)
    if  1-abs(Ci(i)/Ci(1)) >= 0.01 %is this the best way to do it?
        tvmax = t(i-1);
        dv = 30 * z1*(1/t(i+1)-1/t(i));
        Pvmax = 1-abs(Ci(i-1)/Ci(1));
        break
    end
    
end

for j = length(Ci):-1:i
    if  abs(Ci(j)/Ci(1)) >= 0.01
        tvmin = t(j);
        Pvmin = abs(Ci(j));
        break
    end
end
 
vmax= z0/tvmax ;
vmin = z0/tvmin;
v = vmax:dv:vmin;
ti = z0./v;

Ci_i = interp1(t,Ci,ti);

dddtCint = my_2FD_non_uniform(ti,Ci_i);
%dddtCint = movmean(dddtCint,3);

[Pest, Pfit] = P_RHS_solve(z0,z1,ti,dddtCint,Pvmax,Pvmin);
P = Pest;

Nv = length(P);
ci0 = c0 * P;
ci = @(ci0,vs,z,t)(ci0*(1-heaviside(z+vs*t)));
ERR = 0;

%%
if Disp == 1
figure('units','normalized','outerposition',[0 0 1 1])
end

stp = fix(length(t)/50);
r = rem(length(t),stp);

for k = 1:stp:length(t)-r
    
    C = 0;
    for i = 1:Nv
        C = ci(ci0(i),v(i),z,t(k)) + C;
    end
    
    err = sum((C-movmean(C_t_z(:,k),length(z)*0.05)').^2); % OVERALL ERROR
      
    ERR = ERR + err;
       
    if Disp == 1
       [~, boundD] = min(abs(C_t_z(:,k)-c0*0.2)) ; [~, boundU] = min(abs(C_t_z(:,k)-c0*0.8));
       %smallErr = sum((C(1:boundD)-movmean(f_C_t_z(1:boundD,k),length(f_z(1:boundD))*0.05)').^2);
       %bigErr = sum((C(boundU:end)-movmean(f_C_t_z(boundU:end,k),length(f_z(boundU:end))*0.05)').^2);
       
       sb1 = subplot(3,3,[1 4 7]) ;
       
       plot(C,z,'--')
       hold on
       title('Concentration profiles')
       xlim([0 c0*1.1]); ylim([min(z),max(z)])
       plot(movmean(C_t_z(:,k),length(z)*0.05)',z)
       scatter([C_t_z(boundD,k) C_t_z(boundU,k)],[z(boundD) z(boundU)])
       legend('reconstructed data','input data','Measure points')
       set(gca,'FontSize',18)
       
       sb2 = subplot(3,3,[2 5 8]) ;
       
       plot(C_v(:,k)/max(C_v(:,k))*vmax,z,'--')
       title('Ks')
       ylim([min(z),max(z)])
       legend('v(C)')
       set(gca,'FontSize',18)
       
       sb3 = subplot(3,3,3);
       scatter(t(k),err)
       xlim([0 t(end)])
       title('mean square difference for all particles')
       legend(['time : ' sprintf('%.0f',t(k)/t(end)*100) ' [%]'])
       set(gca,'FontSize',18)
       hold on

       if k ==1
           
       sb4 = subplot(3,3,6);
       plot(t,Ci,'--')
       hold on
       plot(ti,Ci_i./max(Ci_i),'LineWidth',2)
       xlim([0 t(end)])
       title('Integration Ci')
       legend('Original','Retained')
       set(gca,'FontSize',18)
       
       sb5 = subplot(3,3,9);
       scatter(v,Pest./max(Pest))
       hold on
       plot(v,Pfit./max(Pfit),'LineWidth', 2)
       title('Particle Velocity Distribution')
       set(gca, 'XLim', [0, get(gca, 'XLim') * [0; 1]])
       ylim([0 1])
       set(gca,'FontSize',18)
       legend('scatter','P fitted')

       end
       
       pause (0.1)
      end
      if err <= 0.001
        break
      end
            
end




end

