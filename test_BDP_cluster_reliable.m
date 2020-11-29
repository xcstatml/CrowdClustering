% BDP_cluster_reliable in a batch mode
clear all;
warning('off','all')

%% setup parameters and generate data
N=12;
K=3;    % number of clusters
GN=N/K;
T=20;   % budget
batch=10;  % batch size 
a0=1;   % prior
b0=1;
at=4;   % true
bt=1;
REP=10;

for rep=1:REP
    C=betarnd(bt,at,N,N);
    clu=[];
    for k=1:K
        clu=[clu;k*ones(GN,1)];
        C(((k-1)*GN+1):(k*GN),((k-1)*GN+1):(k*GN))=betarnd(at,bt,GN,GN);
    end
    C=triu(C,1)+triu(C,1)';
    C(logical(eye(size(C))))=1; % C: matrix of similarity parameters
    
    ndx_okg = BDP_cluster_reliable(a0, b0, T, K, C, batch); % clusters
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