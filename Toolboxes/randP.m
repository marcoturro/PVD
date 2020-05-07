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
close all
x01 = linspace(0,1,refinement);
xpdf = linspace(-35,15,refinement);
muev = 1; bev = 2;
p = evpdf(xpdf,muev,bev);
p = p/sum(p);
ran = zeros(1,50); vals = [-28 -30 -19 -2 -12 -16 -1 0 -3 -10 1 4 9 5 8];
for i = 1:50
ran(i) = vals(randi(numel(vals)));
end

%figure('units','normalized','outerposition',[0 0 1 1])

for r = 2:randi([5 10])
    munorm(r) = ran(randi(numel(ran)));
    if munorm(r) == munorm(r-1)
        munorm(r) = ran(randi(numel(ran)));
    end
    if munorm(r) == munorm(r-1)
        munorm(r) = ran(randi(numel(ran)));
    end
    if munorm(r) < -8
        signorm = (3-1.5)*rand(1) + 1;
        ptmp =  normpdf(xpdf,muev+munorm(r),signorm)/(randi([2 5]));
    else
        signorm = (1-0.5)*rand(1) + 1;
        ptmp = normpdf(xpdf,muev+munorm(r),signorm);
    end
    p = p + ptmp;
end
p = movmean(p,ceil(refinement/30));
%p = interp1(xpdf,p,x01)
p = p/sum(p);

plot(xpdf,p)
end