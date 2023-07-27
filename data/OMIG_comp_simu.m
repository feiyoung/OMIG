% Use a setting to compare OMIG, OrMIG and OfMIG in terms of imputation
% accuracy and computational time/memory usage.
profile on -memory
program.m
profile viewer
profile report


addpath(genpath('./OMIGtool')); ## OMIGtool toolbox can be downloaded at https://github.com/feiyoung/OMIG/tree/main/OMIG_Matlab/OMIGtool
clear
% Given p, test n
n_set = [500, 1000, 3000, 10000, 30000, 50000];
% n_set = [400, 800, 1000];
p=1000;q= 4;type1 = 'norm_pois';
mis_vec = [0.1, 0.2, 0.3] + 0.2; rho = 8;

Nn = length(n_set);
N_rep = 5;
timeArray = zeros(3, Nn,N_rep); % 3 methods, Nn grids of n; N_rep: repeat N_rep times.
aeArray = zeros(3, 2, Nn, N_rep);
naeArray = zeros(3, 2, Nn, N_rep);
maxIter = 10;
for i= 1:Nn 
    % i = 1;
    fprintf('i = %d \n', i);
    n = n_set(i);
    nbatch1 = floor(n*0.5); nbatch = floor((n-nbatch1)/5); nleast = 50;
    idxList = batchGroup(n, nbatch1, nbatch, nleast);
    n_id = length(idxList);
    for ii = 1:N_rep
        % ii = 1;
         fprintf('i = %d, ii= %d \n', i, ii);
         i_rep = (ii-1)*Nn + ii;
        [Xmis, X, H, Bm, group]  =  gendata_mis(i_rep, n, p, type1, q, mis_vec,  [rho, 0.1], 'MAR');

        type = cell(2,2);
        type{1,1} = 'normal'; type{1,2} = 'identity';
        type{2,1} = 'poisson'; type{2,2} = 'log';
        
        % OrMIG
        tic; 
        [hX, hD, hHm, hBm,cvVals, history] = OrMIG(Xmis, group, type, q, 1e-20, maxIter, 1, 1, 1);
        time = toc;
        timeArray(1,i, ii) = time;
        aeArray(1,:,i,ii) = AE(hX, X, Xmis, group);
        naeArray(1,:,i,ii) =  NAE(hX, X, Xmis, group);
        
        % OfMIG
        hX_of = Xmis;
        tic;
        for idx=1:n_id
             [hX_of(idxList{idx},:), hD, hHm, hBm,cvVals, history] = OrMIG(Xmis(idxList{idx},:), group, type, q, 1e-20, maxIter, 1, 1, 1);
        end
        time2 = toc;
        timeArray(2,i, ii) = time2;
        aeArray(2,:,i,ii) = AE(hX_of, X, Xmis, group);
        naeArray(2,:,i,ii) =  NAE(hX_of, X, Xmis, group);
        
        % OMIG
        tic;
        lambda = 0; verbose = 1; parallel = 1; % Bm0 = Bm; H0 = H;
        [X_imp, hXall, hD, Hms, Bms, NAE_mat, corB_vec, corH_vec, errD] = OMIG(Xmis, group, type, q, ...
    idxList, lambda,verbose,parallel, [], [], []);
        time3 = toc;
        timeArray(3,i, ii) = time3;
        aeArray(3,:,i,ii) = AE(X_imp, X, Xmis, group);
        naeArray(3,:,i,ii) =  NAE(X_imp, X, Xmis, group);
        
    end
    
    save ./simu_R1/simuData_Mat/Mdata_R1/Mdata_compOMIG/OMIG_3version_comp.mat timeArray aeArray naeArray
end 
save ./simu_R1/simuData_Mat/Mdata_R1/Mdata_compOMIG/OMIG_3versionp1000_compV2.mat timeArray aeArray naeArray


% Given n, test p

p_set = [500, 1000, 3000, 10000, 30000, 50000];
% p_set = [400, 800, 1000];
n=1000;
Nn = length(p_set);
q= 4;type1 = 'norm_pois';
mis_vec = [0.1, 0.2, 0.3] + 0.2; rho = 8;

nbatch1 = floor(n*0.5); nbatch = floor((n-nbatch1)/5); nleast = 50;
idxList = batchGroup(n, nbatch1, nbatch, nleast);
n_id = length(idxList);
    
Nn = length(n_set);
N_rep = 5;
timeArray = zeros(3, Nn,N_rep); % 3 methods, Nn grids of n; N_rep: repeat N_rep times.
aeArray = zeros(3, 2, Nn, N_rep);
naeArray = zeros(3, 2, Nn, N_rep);
maxIter = 10;
for i= 1:Nn 
    % i = 1;
    fprintf('i = %d \n', i);
    p = p_set(i);
    for ii = 1:N_rep
        % ii = 1;
         fprintf('i = %d, ii= %d \n', i, ii);
        i_rep = (ii-1)*Nn + ii;
        [Xmis, X, H, Bm, group]  =  gendata_mis(i_rep, n, p, type1, q, mis_vec,  [rho, 0.1], 'MAR');

        type = cell(2,2);
        type{1,1} = 'normal'; type{1,2} = 'identity';
        type{2,1} = 'poisson'; type{2,2} = 'log';
        
        % OrMIG
        tic; 
        [hX, hD, hHm, hBm,cvVals, history] = OrMIG(Xmis, group, type, q, 1e-20, maxIter, 1, 1, 1);
        time = toc;
        timeArray(1,i, ii) = time;
        aeArray(1,:,i,ii) = AE(hX, X, Xmis, group);
        naeArray(1,:,i,ii) =  NAE(hX, X, Xmis, group);
        
        % OfMIG
        hX_of = Xmis;
        tic;
        for idx=1:n_id
             [hX_of(idxList{idx},:), hD, hHm, hBm,cvVals, history] = OrMIG(Xmis(idxList{idx},:), group, type, q, 1e-20, maxIter, 1, 1, 1);
        end
        time2 = toc;
        timeArray(2,i, ii) = time2;
        aeArray(2,:,i,ii) = AE(hX_of, X, Xmis, group);
        naeArray(2,:,i,ii) =  NAE(hX_of, X, Xmis, group);
        
        % OMIG
        tic;
        lambda = 0; verbose = 1; parallel = 1; % Bm0 = Bm; H0 = H;
        [X_imp, hXall, hD, Hms, Bms, NAE_mat, corB_vec, corH_vec, errD] = OMIG(Xmis, group, type, q, ...
    idxList, lambda,verbose,parallel, [], [], []);
        time3 = toc;
        timeArray(3,i, ii) = time3;
        aeArray(3,:,i,ii) = AE(X_imp, X, Xmis, group);
        naeArray(3,:,i,ii) =  NAE(X_imp, X, Xmis, group);
        
    end
    
    save ./simu_R1/simuData_Mat/Mdata_R1/Mdata_compOMIG/OMIG_3version_n1000comp.mat timeArray aeArray naeArray
end 
save ./simu_R1/simuData_Mat/Mdata_R1/Mdata_compOMIG/OMIG_3version_n1000comp.mat timeArray aeArray naeArray

