clear;
%cd /home/ligz/LiuWei_files/OMIG
%cd /home/liuwei/Rfiles/OMIG
addpath(genpath('./OMIGtool'));
load DLPFC.mat
Xmis = Dat;
max(max(Dat))

id_var = find(std(Xmis, 'omitnan') ~= 0);
Xmis = Xmis(:, id_var); X = Dat(:, id_var);
[n, p] = size(Xmis);
% generate missing id for data
% 20% missing data
mis_rate_vec = [0.1, 0.2, 0.3, 0.4, 0.5];
nmis = length(mis_rate_vec);
idna10_cell = cell(nmis, 1);
N = 10;
for imis = 1: nmis
    mis_rate = mis_rate_vec(imis);
    n_mis = (n*p)*mis_rate;
    idna10Mat = zeros(n_mis, N);
    for i = 1:N
    rng(i);
    id_na = datasample(1:(n*p), n_mis);
    idna10Mat(:, i) = id_na;
    end
    idna10_cell{imis} = idna10Mat;
end
idna10Mat(1:5, 1:4)
length(id_na) / (n*p)
Xmis(id_na) = NaN;

group = group(id_var);
type = cell(2,2);
type{1,1} = 'normal'; type{1,2} = 'identity';
type{2,1} = 'poisson'; type{2,2} = 'log';
save ./RealData/genBrainMatDat/misRateAllData X idna10_cell group type;



%% Run OrMIG
load misRateAllData.mat
q= 4; N = 10;
n_mis = length(idna10_cell);
type = cell(2,2);
type{1,1} = 'normal'; type{1,2} = 'identity';
type{2,1} = 'poisson'; type{2,2} = 'log';
verbose = 1; lambda =1
timeVec = zeros(N,n_mis);
naeMat = zeros(N,2, n_mis);
aeMat = naeMat;
for imis = 1:n_mis
    % imis = 2
    fprintf('imis = %d \n', imis);
    idna10Mat = idna10_cell{imis};
    size(idna10Mat)
    
    for i = 1:N
        % i = 1;
        length(idna10Mat(:, i)) / (n*p)
        fprintf('i = %d \n', i);
        Xmis = X;
        Xmis(idna10Mat(:, i)) = NaN;
        tic; 
        [hX, hD, hHm, hBm,cvVals, history] = OrMIG(Xmis, group, type, q, 1e-4, 20, lambda, verbose, 1);
        time = toc
        timeVec(i, imis) = time;
        naeMat(i,:, imis) = NAE(hX, X, Xmis, group);
        aeMat(i,:, imis) = AE(hX, X, Xmis, group);
    end
end

save ./RealData/resBrainDat/time_naeMat_gfmImpute_brain76_miss5 naeMat aeMat timeVec;


%% Run the method PLAIS
cd /home/ligz/LiuWei_files/OMIG
% cd /home/liuwei/Rfiles/OMIG
addpath(genpath('./collectivemc'));

load misRateAllData.mat

qVec = zeros(N, n_mis);
naeMat = zeros(N,2, n_mis);
aeMat = naeMat;
for imis = 1:n_mis
    % imis = 2
    fprintf('imis = %d \n', imis);
    idna10Mat = idna10_cell{imis};
    size(idna10Mat)
    
    for i = 1:N
        % i = 1;
        length(idna10Mat(:, i)) / (n*p)
        fprintf('i = %d \n', i);
        Xmis = X;
        Xmis(idna10Mat(:, i)) = NaN;
        Xmis = double(Xmis);
        [n, p] = size(Xmis);
        q = 5;


        % datatypes = {{'gaussian', [.5,  1.]}, {'poisson', 4}, {'bernoulli', 0.5}};
        datatypes = {{'gaussian', [.5,  1.]}, {'poisson', 4}}

        nb_cols = [sum(group==1), sum(group==2)]; % p; %

        ranks = q;

        source1 =Xmis;
        source1(isnan(source1)) = 0;
        [m, n] = size(source1);
        maxS = max(source1(:));
        source_gaussian1 = source1 / maxS;
        source_gaussian1(isnan(source1)) = 0;
        mind = isnan(Xmis);
        mask1 = ~(mind==1);

        [row1, col1, val1] = find(source_gaussian1);

        % cold-start for Poisson data
        [row_gaussian, col_gaussian, data_source_gaussian] = find(source_gaussian1);
        [m_gauss n_gauss] = size(source_gaussian1);
        source_gaussian = source_gaussian1;

        % lambda
        [~, grad_source1] = LikelihoodAndGradLikelihood(source_gaussian, source_gaussian, datatypes, nb_cols, mask1);
        [~, lambda, ~] = lansvd(grad_source1, 1);

        % optimization

        % fprintf("---------------------------------------------------------------\n");
        % fprintf("------------------- PLAIS-IMPUTE on collective cold ----------------\n");
        % fprintf("---------------------------------------------------------------\n\n\n");

        para1.maxIter = 200;
        para1.tol = 1e-5;
        para1.decay = 0.2;
        para1.exact = 0;
        para1.maxR = 10;
        %lambda = 1e-5;

        [U1, S1, V1, output1 ] = PLAISImpute(source_gaussian, lambda, para1, datatypes, nb_cols, mask1);
        qVec(i, imis) = size(U1,2);
        sln_AISoftImpute1 = U1 * S1* V1';
        norm(X(mind) - sln_AISoftImpute1(mind)*maxS)^2/ sum(X(mind).*X(mind))
        hX = sln_AISoftImpute1 * maxS;

       naeMat(i,:,imis) = NAE(hX, X, Xmis, group);
       aeMat(i,:,imis) = AE(hX, X, Xmis, group);
    end
end
save ./RealData/resBrainDat/PLAIS_AE_naeMat naeMat aeMat qVec