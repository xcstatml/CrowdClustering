% BDP_cluster_unreliable in a batch mode
clear all;
warning('off','all')


%% setup parameters and generate data
N=12;
K=3;	% number of clusters
GN=N/K;
W=3;    % number of workers
T=20;	 % budget
batch=100;    % batch size 
a0=1;   % prior
b0=1;
c0=4;
d0=1;
at=4;   % true
bt=1;
ct=8;
dt=1;
REP=10;

for rep=1:REP
    CR=betarnd(bt,at,N,N);
    clu=[];
    for k=1:K
        clu=[clu;k*ones(GN,1)];
        CR(((k-1)*GN+1):(k*GN),((k-1)*GN+1):(k*GN))=betarnd(at,bt,GN,GN);
    end
    CR=triu(CR,1)+triu(CR,1)';
    CR(logical(eye(size(CR))))=1;   % CR: matrix of similarity parameters
    
    for w=1:W
        WR(w)=0.65+0.05*w;  % WR: workers' reliability
    end
    
    ndx_okg = BDP_cluster_unreliable(a0, b0, c0, d0, T, K, CR, WR, batch);  % clusters
    for t=1:T        
        IN=perms(1:K);
        for s=1:size(IN,1)
            may=reshape(ones(GN,1)*IN(s,:),1,GN*K);
            may_okg(s)=length(find(ndx_okg(t,:)-may==0));
        end
        acc_okgr(t,rep)=max(may_okg)/N;
    end
end
acc=mean(acc_okgr,2)    % clustering accuracy