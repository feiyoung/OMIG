clear;
cd /home/liuwei/Rfiles/gfmImpute
addpath(genpath('./matFiles'));
addpath('./glmnet');

%% Generate continuous data with only one variable type
clear;
i = 4; p = 100; n = 100; q = 6;type = 'heternorm';
[Xmis, X, H, Bm, group] = gendata_mis(i, n, p, type, q, [0.1, 0.2, 0.3], 8, 'MNAR');

% unified function
group = [ones(1,p)]; % full-one vector indicates all variables belong to the same type.
type = cell(1,2);
type{1,1} = 'normal';  % the type is 'normal'
type{1,2} = 'identity'; % the link funciton is 'identity'
q= 6; lambda=0; verbose=1;
% imputation using gfmImpute
tic; 
[hX, hD, hHm, hBm,cvVals, history] = OrMIG(Xmis, group, type, q, 1e-6, 10, lambda, verbose, 0);
time = toc
% [measurefun([ones(n,1) H], hHm), ...
% measurefun(Bm, hBm)]
NAE(hX, X, Xmis, group)

%% Generate mixed-type  data with  two variable types including count and binary variables.
i = 1; p = 200; n = 200; q=4; vtype = 'pois_bino';
[Xmis, X, H, Bm, group] = gendata_mis(i, n, p, vtype, q, [0.1, 0.2, 0.3], 8, 'MAR');

% unified function
type = cell(2,2);
type{1,1} = 'poisson'; type{1,2} = 'log';
type{2,1} = 'binomial';  type{2,2} = 'probit'; 
q= 4; lambda=0; verbose=1;
% imputation using gfmImpute
tic; 
[hX, hD, hHm, hBm,cvVals, history] = OrMIG(Xmis, group, type, q, 1e-6, 10, lambda, verbose, 0);
time = toc
[measurefun([ones(n,1) H], hHm), ...
measurefun(Bm, hBm)]
NAE(hX, X, Xmis, group)

%% Generate mixed-type data with  two variable types including continuous and binary variables.
i = 1; p = 200; n = 200; q=4; type = 'norm_bino';
[Xmis, X, H, Bm, group] = gendata_mis(i, n, p, type, q, 0.4,  8, 'MCAR');

% unified function
type = cell(2,2);
type{1,1} = 'normal'; type{1,2} = 'identity';
type{2,1} = 'binomial';  type{2,2} = 'probit'; 
q= 6; lambda=0.1; verbose=1;
% imputation using gfmImpute
tic; 
[hX, hD, hHm, hBm,cvVals, history] = OrMIG(Xmis, group, type, q, 1e-6, 10, lambda, verbose, 0);
time = toc
[measurefun([ones(n,1) H], hHm), ...
measurefun(Bm, hBm)]
NAE(hX, X, Xmis, group)


%% Generate mixed-type  data with  three variable types including continuous, count and binary variables.
n = 200; p=200; i = 4;q=3; type = 'npb';
[Xmis, X, H, Bm, group]  = gendata_mis(i, n, p, type, q, 0.3, 6);
max(max(X))
% unified functions test
id_var = find(std(Xmis, 'omitnan') ~= 0);
Xmis = Xmis(:, id_var); X = X(:, id_var);
group = group(id_var);
type = cell(3,2);
type{1,1} = 'normal'; type{1,2} = 'identity';
type{2,1} = 'poisson'; type{2,2} = 'log';
type{3,1} = 'binomial';  type{3,2} = 'logit';
q= 3;
verbose = 1; lambda =10;
% imputation using gfmImpute
tic; 
[hX, hD, hHm, hBm,cvVals, history] = OrMIG(Xmis, group, type, q, 1e-6, 10, lambda, verbose, 0);
time = toc
% [measurefun([ones(n,1) H], hHm), ...
% measurefun(Bm, hBm)]
NAE(hX, X, Xmis, group)
% Select 
q_set = 2:4;
lambda = 0;
parallel = 1; verbose = 1;
[icMat] = selectFacNumber(Xmis,group, type, q_set,lambda, parallel, verbose);
icMat
[~, index_min] = min(icMat(:,1)); 
fprintf('----------The best number of factors by IC criteria is: %d -----\n',q_set(index_min) )
[~, index_min_pc] = min(icMat(:,2)); 
fprintf('----------The best number of factors by PC criteria is: %d -----\n',q_set(index_min_pc) )


