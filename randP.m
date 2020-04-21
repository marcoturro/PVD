function [x01, p] = randP(refinement)
%%
% Input the refinement of the PVD and output a PVD made out of 3 to 4
% peeks defined by the random term "r" in the loop and by an initial GUMMEL
% or EXTREME VALUE distribution to create a left skewness in the
% Distribution.
%
% The term used for the initial EV distribution are hand tuned to be
% similar to known PSD taken from GILLARD () and from known sediment
% analysis of CCFZ samples
%

x01 = linspace(0,1,refinement);
xpdf = linspace(-16,5,refinement);
muev = -2; bev = 2;
p = evpdf(xpdf,muev,bev);
p = p/sum(p);
ran = zeros(1,50); vals = [-5 -4 -3 -2 -1 0 -3 -1 1 3 4 3 4 5];
for i = 1:50
ran(i) = vals(randi(numel(vals)));
end

figure('units','normalized','outerposition',[0 0 1 1])

for r = 2:randi([4 6])
    signorm = (1.5-1)*rand(1) + 1;
    munorm(r) = ran(randi(numel(ran)));
    if munorm(r) == munorm(r-1)
        munorm(r) = ran(randi(numel(ran)));
    end
    if munorm(r) == munorm(r-1)
        munorm(r) = ran(randi(numel(ran)));
    end
    ptmp = normpdf(xpdf,muev+munorm(r),signorm);
    beta = (1.2-0.4).*rand(1) + 0.4;
    p = p + ptmp/(beta*sum(ptmp));
    plt = plot(x01,p,'LineWidth',3);
%     pause(0.1)
%     hold on
%     legend()
end

p(1)=0;
p = p/sum(p);
end