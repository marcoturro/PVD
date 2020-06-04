function NoiseAnalysis(data,NPC)
 
figure

M = 10;    % window length = embedding dimension
t = data.t;
if length(data.C) == t
    data.C = data.C';
end
N = length(data.C(2,:));   % length of generated time series
T = data.t(2)-data.t(1);    % period length of sine function
% stdnoise = 1; % noise-to-signal ratio

X = normalize(data.C(floor(length(data.C))/2,:));

covX = xcorr(X,M-1,'unbiased');
Ctoep=toeplitz(covX(M:end));

figure(1);
set(gcf,'name','Covariance matrix');
clf;
imagesc(Ctoep);
axis square
colorbar

Y=zeros(N-M+1,M);
for m=1:M
  Y(:,m) = X((1:N-M+1)+m-1);
end;
Cemb=Y'*Y / (N-M+1);

% figure(2);
% set(gcf,'name','Covariance matrix');
% clf;
% imagesc(Cemb);
% axis square
% colorbar

%%%%% ------- %%%%%


C = Ctoep;
% C = Cemb;

[RHO,LAMBDA] = eig(C);
LAMBDA = diag(LAMBDA);               % extract the diagonal elements
[LAMBDA,ind]=sort(LAMBDA,'descend'); % sort eigenvalues
RHO = RHO(:,ind);                    % and eigenvectors
% 
% figure(3);
% set(gcf,'name','Eigenvectors RHO and eigenvalues LAMBDA')
% clf;
% subplot(3,1,1);
% plot(LAMBDA,'o-');
% subplot(3,1,2);
% plot(RHO(:,1:2), '-');
% legend('1', '2');
% subplot(3,1,3);
% plot(RHO(:,3:4), '-');
% legend('3', '4');

%%%%% ------- %%%%%

PC = Y*RHO;

figure(4);
set(gcf,'name','Principal components PCs')
clf;
for m=1:4
  subplot(4,1,m);
  plot(t(1:N-M+1),PC(:,m),'k-');
  ylabel(sprintf('PC %d',m));
  ylim([-10 10]);
end;


RC=zeros(N,M);
for m=1:M
  buf=PC(:,m)*RHO(:,m)'; % invert projection
  buf=buf(end:-1:1,:);
  for n=1:N % anti-diagonal averaging
    RC(n,m)=mean( diag(buf,-(N-M+1)+n) );
  end
end;

figure(5);
set(gcf,'name','Reconstructed components RCs')
clf;
for m=1:4
  subplot(4,1,m);
  plot(t,RC(:,m),'r-');
  ylabel(sprintf('RC %d',m));
  ylim([-1 1]);
end;

%%%%% ------- %%%%%
%%

figure(6);
set(gcf,'name','Original time series X and reconstruction RC')
clf;
subplot(NPC,1,1)
plot(t,X,'b-');
legend('Original');


x = zeros(NPC,length(t));

for i = 1:NPC% nb of PC
    x(i,:) = sum(RC(:,i),2); % 3th and 4th component
end

Fs = 1/T;
for i = 1:NPC
subplot(5,1,i)
plot(t,x(i,:));

end

%%%%% ------- %%%%%

figure
for i = 1:NPC-1
subplot(2,2,i)
psdest = psd(spectrum.periodogram,x(i+1,:),'Fs',Fs,'NFFT',length(x));
[~,I] = max(psdest.Data);
plot(psdest)
end

end