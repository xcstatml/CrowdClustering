function ndx_okg = BDP_cluster_reliable(a0, b0, T, K, C, batch)

% (a0,b0): parameters of prior Beta distributions
% T: budget
% K: number of clusters
% C: matrix of similarity parameters
% batch: batch size

N=size(C,1);
GN=N/K;

A_okg=a0*ones(N,N);
B_okg=b0*ones(N,N);
Wokg=A_okg./(A_okg+B_okg);
Wokg(logical(eye(size(Wokg))))=1;

for t=1:T
    OKG=1000*ones(N);
    COMOKG=[];
    for i=1:(N-1)  
        Vokg=zeros(1, N-i);
        for j=(i+1):N
            A=A_okg;
            B=B_okg;
            W1=Wokg;
            W2=Wokg;
            W1(i,j)=(A(i,j)+1)/(A(i,j)+B(i,j)+1);
            W1(j,i)=W1(i,j);
            W2(i,j)=A(i,j)/(A(i,j)+B(i,j)+1);
            W2(j,i)=W2(i,j);
            [ndx1,Pi1,cost1]= grPartition(W1,K,30);
            [ndx2,Pi2,cost2]= grPartition(W2,K,30);
            Vokg(j-i)=min(cost1,cost2);
        end
        COMOKG=[COMOKG,Vokg];
    end

    for i=1:(N-1)
        st=(2*N-i)*(i-1)/2+1;
        en=(2*N-i-1)*i/2;
        OKG(i,(i+1):N)=COMOKG(st:en);
    end

    [AA_OKG,I_OKG]=sort(OKG(:));
    [x2,y2]=ind2sub(size(OKG),I_OKG);

    for s=1:batch
        l=binornd(1,C(x2(s),y2(s)));    % similarity label
        A_okg(x2(s),y2(s))=A_okg(x2(s),y2(s))+l;
        A_okg(y2(s),x2(s))=A_okg(x2(s),y2(s));
        B_okg(x2(s),y2(s))=B_okg(x2(s),y2(s))+1-l;
        B_okg(y2(s),x2(s))=B_okg(x2(s),y2(s));
    end

    Wokg=A_okg./(A_okg+B_okg);
    Wokg(logical(eye(size(Wokg))))=1;
    [ndx_okg(t,:),Pi,cost_okg(t)]= grPartition(Wokg,K,30);
end