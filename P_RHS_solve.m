function [P] = P_RHS_solve(z0,z1,ti,dddtCint,Pvmax,Pvmin)
%%
v = z0./ti;
a = z1/z0;
v_star = v*a;
A = zeros(length(ti));

rpgddt = fitrgp(ti',dddtCint','Basis','linear',...
      'FitMethod','exact','PredictMethod','exact'); %,'KernelFunction','exponential'); for rougher datasets?
%dddtCint = resubPredict(rpgddt)';

for i = 1:length(v)
    
    idx = find(v<v_star(i), 1 ); %this creates issues
    coef1 = abs(v_star(i) - v(idx));
    coef2 = abs(v_star(i) - v(idx-1));
    dv = abs(v(idx)- v(idx-1));
    if(isempty(dv))
        coef1 = 0; %coef1/dv;
        coef2 = 0;%coef2/dv;
    else
        coef1 =  coef1/dv;
        coef2 = coef2/dv;
    end
    
    A(i,idx) = - coef2 * a^2;
    A(i,idx-1) = - coef1 * a^2;
    A(i,i) = 1;
    
end
A(1,:) = 0;
A(1,1) = 1;
A(end,:)=0;
A(end,end)=1;

rhs = (dddtCint./v.^3)';
rhs(1) = Pvmax ; %vmax condition
rhs(end) = Pvmin ; %vmin condition
P = linsolve(A,rhs);

%interp1([0 1],[0 vs(1)],0:dv:vs(1))
%P(end+1) = P(end)/v(end);
%v(end+1) = 0;


%P = flip(P');
%v = flip(v');
%dv = mean(diff(v));

% v(P<0) = [];
% P(P<0) = [];

% Pi = interp1([0 v(1)],[0 P(1)],0:dv:v(1));
% P = [Pi P(2:end)];
% vi = interp1([0 v(1)],[0 v(1)],0:dv:v(1));
% v = [vi v(2:end)];

%%
%rpgP = fitrgp(v',P','Basis','linear',...
 %     'FitMethod','exact','PredictMethod','exact'); %,'KernelFunction','exponential'); for rougher datasets?

%Pfit = resubPredict(rpgP);
%Pfit = P;
%Pfit = Pfit./sum(Pfit);
%P = P./sum(P);

end
