function ndx_okg = BDP_cluster_unreliable(a0, b0, c0, d0, T, K, CR, WR, batch)

% (a0,b0),(c0,d0): parameters of prior Beta distributions
% T: budget
% K: number of clusters
% CR: matrix of similarity parameters
% WR: workers' reliability
% batch: batch size

N=size(CR,1);
GN=N/K;
W=size(WR,2);

A_okg=a0*ones(N,N);
B_okg=b0*ones(N,N);
C_okg=c0*ones(W,1);
D_okg=d0*ones(W,1);
Wokg=A_okg./(A_okg+B_okg);
Wokg(logical(eye(size(Wokg))))=1;

for t=1:T        
    for w=1:W
        OKG(:,:,w)=1000*eye(N);
        COMOKG=[];
        for i=1:(N-1)    % parallel
            Vokg=zeros(1, N-i);
            for j=(i+1):N
                A=A_okg;
                B=B_okg;
                C=C_okg;
                D=D_okg;
                W1=Wokg;
                W2=Wokg;
                E1t=A(i,j)*((A(i,j)+1)*C(w)+B(i,j)*D(w))/(A(i,j)+B(i,j)+1)/(A(i,j)*C(w)+B(i,j)*D(w));
                E1t2=A(i,j)*(A(i,j)+1)*((A(i,j)+2)*C(w)+B(i,j)*D(w))/(A(i,j)+B(i,j)+1)/(A(i,j)+B(i,j)+2)/(A(i,j)*C(w)+B(i,j)*D(w));
                E0t=A(i,j)*(B(i,j)*C(w)+(A(i,j)+1)*D(w))/(A(i,j)+B(i,j)+1)/(B(i,j)*C(w)+A(i,j)*D(w));
                E0t2=A(i,j)*(A(i,j)+1)*(B(i,j)*C(w)+(A(i,j)+2)*D(w))/(A(i,j)+B(i,j)+1)/(A(i,j)+B(i,j)+2)/(B(i,j)*C(w)+A(i,j)*D(w));
                aa=E1t*(E1t-E1t2)/(E1t2-E1t^2);
                bb=(1-E1t)*(E1t-E1t2)/(E1t2-E1t^2);
                W1(i,j)=aa/(aa+bb);
                W1(j,i)=W1(i,j);
                aa=E0t*(E0t-E0t2)/(E0t2-E0t^2);
                bb=(1-E0t)*(E0t-E0t2)/(E0t2-E0t^2);
                W2(i,j)=aa/(aa+bb);
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
            OKG(i,(i+1):N,w)=COMOKG(st:en);
        end
    end
    
    [AA_OKG,I_OKG]=sort(OKG(:));
    [x2,y2,z2]=ind2sub(size(OKG),I_OKG);

    for s=1:batch
        i=x2(s);
        j=y2(s);
        w=z2(s);
        p=WR(z2(s))*CR(x2(s),y2(s))+(1-WR(z2(s)))*(1-CR(x2(s),y2(s)));
        l=binornd(1,p); % similarity label
        A=A_okg;
        B=B_okg;
        C=C_okg;
        D=D_okg;
        E1t=A(i,j)*((A(i,j)+1)*C(w)+B(i,j)*D(w))/(A(i,j)+B(i,j)+1)/(A(i,j)*C(w)+B(i,j)*D(w));
        E1t2=A(i,j)*(A(i,j)+1)*((A(i,j)+2)*C(w)+B(i,j)*D(w))/(A(i,j)+B(i,j)+1)/(A(i,j)+B(i,j)+2)/(A(i,j)*C(w)+B(i,j)*D(w));
        E0t=A(i,j)*(B(i,j)*C(w)+(A(i,j)+1)*D(w))/(A(i,j)+B(i,j)+1)/(B(i,j)*C(w)+A(i,j)*D(w));
        E0t2=A(i,j)*(A(i,j)+1)*(B(i,j)*C(w)+(A(i,j)+2)*D(w))/(A(i,j)+B(i,j)+1)/(A(i,j)+B(i,j)+2)/(B(i,j)*C(w)+A(i,j)*D(w));
        E1r=C(w)*(A(i,j)*(C(w)+1)+B(i,j)*D(w))/(C(w)+D(w)+1)/(A(i,j)*C(w)+B(i,j)*D(w));            
        E1r2=C(w)*(C(w)+1)*(A(i,j)*(C(w)+2)+B(i,j)*D(w))/(C(w)+D(w)+1)/(C(w)+D(w)+2)/(A(i,j)*C(w)+B(i,j)*D(w));          
        E0r=C(w)*(B(i,j)*(C(w)+1)+A(i,j)*D(w))/(C(w)+D(w)+1)/(B(i,j)*C(w)+A(i,j)*D(w));            
        E0r2=C(w)*(C(w)+1)*(B(i,j)*(C(w)+2)+A(i,j)*D(w))/(C(w)+D(w)+1)/(C(w)+D(w)+2)/(B(i,j)*C(w)+A(i,j)*D(w));          
        if l==1
            A_okg(x2(s),y2(s))=E1t*(E1t-E1t2)/(E1t2-E1t^2);
            A_okg(y2(s),x2(s))=A_okg(x2(s),y2(s));
            B_okg(x2(s),y2(s))=(1-E1t)*(E1t-E1t2)/(E1t2-E1t^2);
            B_okg(y2(s),x2(s))=B_okg(x2(s),y2(s));
            C_okg(z2(s))=E1r*(E1r-E1r2)/(E1r2-E1r^2);
            D_okg(z2(s))=(1-E1r)*(E1r-E1r2)/(E1r2-E1r^2);
        else
            A_okg(x2(s),y2(s))=E0t*(E0t-E0t2)/(E0t2-E0t^2);
            A_okg(y2(s),x2(s))=A_okg(x2(s),y2(s));
            B_okg(x2(s),y2(s))=(1-E0t)*(E0t-E0t2)/(E0t2-E0t^2);
            B_okg(y2(s),x2(s))=B_okg(x2(s),y2(s));
            C_okg(z2(s))=E0r*(E0r-E0r2)/(E0r2-E0r^2);
            D_okg(z2(s))=(1-E0r)*(E0r-E0r2)/(E0r2-E0r^2);
        end
    end

    Wokg=A_okg./(A_okg+B_okg);
    Wokg(logical(eye(size(Wokg))))=1;
    [ndx_okg(t,:),Pi,cost_okg(t)]= grPartition(Wokg,K,30);
end