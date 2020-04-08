function [v, P] = PVD_solve(dat)
%%
clearvars -except dat Disp

f_t = dat.t;
f_z = dat.z;
f_Ci = dat.Ci;
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

% order = 3; frame = 13;
% f_Ci(2:end) = sgolayfilt(f_Ci(2:end),order,frame);
% f_Ci(2:end)= wdenoise(f_Ci(2:end),4);
for i = 1:length(f_Ci)
    if  1-abs(f_Ci(i)/f_Ci(1)) >= 0.00001 %is this the best way to do it?
        tvmax = f_t(i-1);
        %dv = 30 * f_z1*(1/f_t(i+1)-1/f_t(i));
        Pvmax = 1-abs(f_Ci(i-1)/f_Ci(1));
        break
    end
    
end

for j = length(f_Ci):-1:i
    if  abs(f_Ci(j)/f_Ci(1)) >= 0.01
        tvmin = f_t(j);
        Pvmin = abs(f_Ci(j));
        break
    end
end
 
vmax= f_z1/tvmax ;
vmin = f_z0/tvmin;
%v = vmax:dv:vmin;
v = linspace(vmax,vmin,30);
ti = f_z0./v;

tint = ti;

Ci_i = interp1(f_t,f_Ci,ti);

dddtCint = my_2FD_non_uniform(ti,Ci_i);
%dddtCint = movmean(dddtCint,3);
[Pest] = P_RHS_solve(f_z0,f_z1,ti,dddtCint,Pvmax,Pvmin);
P = Pest;

Nv = length(P);
ci0 = f_c0 * P;
ci = @(ci0,vs,z,t)(ci0*(1-heaviside(z+vs*t)));
ERR = 0;

% %%
% if Disp == 1
% figure('units','normalized','outerposition',[0 0 1 1])
% end

% 
% 
% for k = 1:length(ti)
%     
%     C = 0;
%     for i = 1:Nv
%         C = ci(ci0(i),v(i),f_z,ti(k)) + C;
%     end
%     
% %        err = sum((C-movmean(f_C_t_z(:,k),length(f_z)*0.05)').^2); % OVERALL ERROR
% %       
% %        ERR = ERR + err;
% %        
%         if Disp == 1
% %        [~, boundD] = min(abs(f_C_t_z(:,k)-f_c0*0.2)) ; [~, boundU] = min(abs(f_C_t_z(:,k)-f_c0*0.8));
% %        smallErr = sum((C(1:boundD)-movmean(f_C_t_z(1:boundD,k),length(f_z(1:boundD))*0.05)').^2);
% %        bigErr = sum((C(boundU:end)-movmean(f_C_t_z(boundU:end,k),length(f_z(boundU:end))*0.05)').^2);
% %        sb1 = subplot(3,2,[1 3 5]) ;
% %        
% %        plot(C,f_z,'--')
% %        hold on
% %        title('Concentration profiles')
% %        xlim([0 f_c0*1.1]); ylim([min(f_z),max(f_z)])
% %        plot(movmean(f_C_t_z(:,k),length(f_z)*0.05)',f_z)
% %        scatter([f_C_t_z(boundD,k) f_C_t_z(boundU,k)],[f_z(boundD) f_z(boundU)])
% %        legend('reconstructed data','input data','Measure points')
% %        set(gca,'FontSize',18)
% %        
% %        sb2 = subplot(3,2,2);
% %        scatter(f_t(k),err)
% %        xlim([0 f_t(end)])
% %        title('mean square difference for all particles')
% %        legend(['time : ' sprintf('%.0f',f_t(k)*100) ' [%]'])
% %        set(gca,'FontSize',18)
% 
%        if k ==1
%            
%        sb3 = subplot(3,2,4);
%        plot(f_t,f_Ci,'--')
%        hold on
%        plot(ti,Ci_i./max(Ci_i),'LineWidth',2)
%        xlim([0 f_t(end)])
%        title('Integration Ci')
%        legend('Original','Retained')
%        set(gca,'FontSize',18)
%        
%        sb4 = subplot(3,2,6);
%        scatter(v,Pest./max(Pest))
%        hold on
%        plot(v,Pfit./max(Pfit),'LineWidth', 2)
%        title('Particle Velocity Distribution')
%       % ylim([0 1])
%        set(gca,'FontSize',18)
%        legend('scatter','P fitted')
% 
%        end
%        pause (0.1)
%        end
% %        if err <= 0.001
% %         break
% %        end
%             
% end
% 
% 


end

