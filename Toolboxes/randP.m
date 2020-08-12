
function [xpdf, p] = randP(refinement)
%%
% Input the refinement of the PVD and output a PVD made out of 3 to 4
% peeks defined by the random term "r" in the loop and by an initial GUMMEL
% or EXTREME VALUE distribution to create a left skewness in the
% Distribution.
%
% The term used for the initial EV distribution are hand tuned to be
% similar to known PSD taken from GILLARD () and from known sediment
% analysis of CCFZ samples
 
s1 = logspace(0,3,refinement*0.4-1)/10^3*0.005; 
s2 = logspace(0,3,refinement*0.5)/10^3;
[~, id] = min(abs(s2-s1(end)));
xpdf = [0 s1 s2(id+1:end)];
p = zeros(1,length(xpdf));
 
for i = 0:randi([20 50],1)
    mui = xpdf(randi([floor(refinement/2) length(xpdf)],1));   
    sigi = randi([200 700])/100000;
    pi = normpdf(xpdf,mui,sigi);
    p = p + pi;

end

p(1) = 0;
p = p/sum(p);
% figure
% bar(log(xpdf),p,'FaceAlpha',0.3,'EdgeColor','none');
% hold on
 
 
end
