function ERR = PVD_Opt(x,dat,Disp)

f_t = dat{1};
f_z = dat{2};
f_Cimavg = dat{3};
f_P0 = dat{4};
f_c0 = dat{5};
f_z0 = dat{6};
f_z1 = dat{7};
f_C_t_z = dat{8};

 if ~exist('Disp','var')
     % third parameter does not exist, so default it to something
      Disp = 0;
 end
 
Nt = x(1);
tSwitch = x(2);
ti = linspace(f_t(1),tSwitch,Nt);
refinement = 200; % input('specify number of points for data interpolation \nA higher value will create a finer result')
try
    start = ti(end); %input
catch
    start = 5;
    ti(1) = start;
    Nt = 2;
end

DT = (f_t(end)/start)^(1/refinement);

for k = Nt+1:(Nt + refinement)
    ti(k) = ti(k-1)*DT;
end

ti = movmean(ti,Nt);
%%

Ci_i = interp1(f_t,f_Cimavg,ti);

order = 3; frame = 15;
Ci_i = sgolayfilt(Ci_i,order,frame);

dddtCint = my_2FD_non_uniform(ti,Ci_i);

[v, P] = P_RHS_solve(f_z0,f_z1,ti,dddtCint,f_P0);

P(P<=0)=0; % This is to be revisited
P = movmean(P,3);
P = P/sum(P);


%% verification


Nv = size(P);
vr = - ones(length(P),1);
ci0 = f_c0 * P;
ci = @(ci0,vs,z,t)(ci0*(1-heaviside(z+vs*t)));
ERR = 0;
boundD = floor(length(f_z)*0.05); boundU = floor(length(f_z)*0.7);
if Disp == 1
figure('units','normalized','outerposition',[0 0 1 1])
end
for k = 1:length(f_t)
    
    C = 0;
    for i = 1:Nv
        C = ci(ci0(i),v(i),f_z,f_t(k)) + C;
    end
       err = sum((C-movmean(f_C_t_z(:,k),length(f_z)*0.05)').^2); %OVERALL ERROR
       
       %err = sum((C(1:boundD)-movmean(f_C_t_z(1:boundD,k),length(f_z(1:boundD))*0.05)').^2);%SMALL PARTICLES ERROR
       
       ERR = ERR + err;
       
       if Disp == 1
       smallErr = sum((C(1:boundD)-movmean(f_C_t_z(1:boundD,k),length(f_z(1:boundD))*0.05)').^2);
       bigErr = sum((C(boundU:end)-movmean(f_C_t_z(boundU:end,k),length(f_z(boundU:end))*0.05)').^2);
       subplot(3,2,[1 3 5]) 
       plot(C,f_z)
       title('Concentration profiles')
       if k == 1; xl = xlim; yl = ylim; end
       xlim(xl); ylim(yl)
       hold on
       plot(movmean(f_C_t_z(:,k),length(f_z)*0.05)',f_z,'-x')
       legend('reconstructed data','input data','Measure point 1','Measure point 2')
       scatter([C(boundD) C(boundU)],[f_z(boundD) f_z(boundU)])
       hold off
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
       scatter(f_t(k),err)
       xlim([0 f_t(end)])
       title('mean square difference for all particles')
       %legend(['time remaining : ' sprintf('%.0f',f_t(end)-f_t(k)) ' [s]'])
       hold(sb2,'on') ; hold(sb3,'on') ; hold(sb4,'on')
       set(gca,'FontSize',18)
       pause (0.05)
       end
            
end

if Disp == 1
    figure
    plot(v,P)
end
end