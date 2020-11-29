function [ndx,Pi,cost]= grPartition(C,k,nrep);
% Inputs:
% C - n by n edge-weights matrix.
% k - desired number of partitions
% nrep - number of repetion for the clustering algorithm
% Outputs:
% ndx - n-vector with the cluster index for every node
% Pi - Projection matrix
% cost - cost of the partition (sum of broken edges)
%
% By Joao Pedro Hespanha, Copyright 2004
% Make C double stochastic
n=size(C,1);
C=C/(max(sum(C)));
C=C+sparse(1:n,1:n,1-sum(C));
% Spectral partition
options.issym=1; % matrix is symmetric
options.isreal=1; % matrix is real
options.tol=1e-6; % decrease tolerance
options.maxit=500; % increase maximum number of iterations
options.disp=0;
[U,D]=eigs(C,k,'la',options); % only compute 'k' largest eigenvalues/vectors
%[U,D]=eigs(C,k,'LR',options); % only compute 'k' largest eigenvalues/vectors
% Clustering -- requires the Statistics Toolbox
[ndx,zs]=kmeans(U,k,'Distance','cosine','Start','sample','Replicates',nrep);
Pi=sparse(1:length(ndx),ndx,1);
cost=full(sum(sum(C))-trace(Pi'*C*Pi));