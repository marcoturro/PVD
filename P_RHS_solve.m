function [v, P] = P_RHS_solve(z0,z1,ti,dddtCint,P0)
%%
v = z0./ti;
a = z1/z0;
v_star = v*a;
A = zeros(length(ti));

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
A(1,2) = 0;

rhs = (dddtCint./v.^3)';
rhs(1) = 0 ; %Inital value &
rhs(end) = 0 ; %Last value problem
P = linsolve(A,rhs);
P(end+1) = P(end)/v(end);
v(end+1) = 0;

P = P/sum(P);

end
