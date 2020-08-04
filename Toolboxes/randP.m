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
p = zeros(1,refinement);
xpdf = [0 logspace(0,2,refinement-1)/10^2];
for i = 0:randi([100 200],1)
    mui = xpdf(randi([1 refinement],1))
    sigi = rand/10;
    pi = normpdf(xpdf,mui,sigi);
    p = p + pi/sum(pi);
end
p(1) = 0;

end