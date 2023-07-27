cd ./gfmImpute

% collectivemc matlab toolbox is downloaed at https://github.com/mzalaya/collectivemc/
addpath(genpath('./collectivemc'));
addpath(genpath('./matFiles'));


% load ./gfmImpute/Rdata/simu3_misRate2mechMCARXlist.mat
% load ./gfmImpute/Rdata/simu3_misRate2mechMARXlist.mat
XmisCell = struct2cell(XmisList);
XCell =  struct2cell(XList);
clear XmisList XList
n  = length(XmisCell);
meaMat = zeros(n, 3);
qVec = zeros(n,1);
for i = 1:n
    % i= 1
Xmis = XmisCell{i};
X = XCell{i};
X = double(X);
Xmis = double(Xmis);
Xmis(Xmis==-2147483648) = NaN;
mind = isnan(Xmis);
sum(mind)

Xmis = double(Xmis);
[n, p] = size(Xmis);
q = 5;


% datatypes = {{'gaussian', [.5,  1.]}, {'poisson', 4}, {'bernoulli', 0.5}};
% datatypes = {{'poisson', 4}, {'bernoulli', 0.8}}
datatypes = {{'gaussian', [.5,  1.]}, {'poisson', 4}, {'bernoulli', 0.5}};

nb_cols = [sum(group==1), sum(group==2), sum(group==3)]; % p; %

ranks = q;

source1 =Xmis;
source1(isnan(source1)) = 0;
[m, n] = size(source1);
maxS = max(source1(:));
source_gaussian1 = source1 / maxS;
source_gaussian1(isnan(source1)) = 0;

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
para1.maxR = 5* q;
%lambda = 1e-5;

[U1, S1, V1, output1 ] = PLAISImpute(source_gaussian, lambda, para1, datatypes, nb_cols, mask1);
qVec(i) = size(U1,2);
sln_AISoftImpute1 = U1 * S1* V1';
norm(X(mind) - sln_AISoftImpute1(mind)*maxS)^2/ sum(X(mind).*X(mind))
hX = sln_AISoftImpute1 * maxS;

meaMat(i,:) = NAE(hX, X, Xmis, group);
% NAE(sln_AISoftImpute1, X / maxS, Xmis/maxS, group)
end

mean(meaMat)

% save ./gfmImpute/Rdata/simu3_misRate2mechMCAR_PLAIS qVec meaMat
% save ./gfmImpute/Rdata/simu3_misRate2mechMAR_PLAIS qVec meaMat
save ./gfmImpute/Rdata/simu3_misRate2mechMNAR_PLAIS qVec meaMat


%% CSDA R1
cd /home/ligz/LiuWei_files/OMIG
addpath(genpath('./collectivemc'));
clear
% load ./simu_R1/simuData_Mat/Rdata/simu1_misRate1mechMCARXlist.mat
load ./simu_R1/simuData_Mat/Rdata/simu3_misRate2mechMARXlist.mat
XmisCell = struct2cell(XmisList);
XCell =  struct2cell(XList);
clear XmisList XList
n  = length(XmisCell);
meaMat = zeros(n, 3);
qVec = zeros(n,1);
for i = 1:n
    % i= 1
Xmis = XmisCell{i};
X = XCell{i};
X = double(X);
Xmis = double(Xmis);
Xmis(Xmis==-2147483648) = NaN;
mind = isnan(Xmis);
sum(mind)

Xmis = double(Xmis);
[n, p] = size(Xmis);
q = 5;


% datatypes = {{'gaussian', [.5,  1.]}, {'poisson', 4}, {'bernoulli', 0.5}};
% datatypes = {{'poisson', 4}, {'bernoulli', 0.8}}
datatypes = {{'gaussian', [.5,  1.]}, {'poisson', 4}, {'bernoulli', 0.5}};

nb_cols = [sum(group==1), sum(group==2), sum(group==3)]; % p; %

ranks = q;

source1 =Xmis;
source1(isnan(source1)) = 0;
[m, n] = size(source1);
maxS = max(source1(:));
source_gaussian1 = source1 / maxS;
source_gaussian1(isnan(source1)) = 0;

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
para1.maxR = 5* q;
%lambda = 1e-5;

[U1, S1, V1, output1 ] = PLAISImpute(source_gaussian, lambda, para1, datatypes, nb_cols, mask1);
qVec(i) = size(U1,2);
sln_AISoftImpute1 = U1 * S1* V1';
norm(X(mind) - sln_AISoftImpute1(mind)*maxS)^2/ sum(X(mind).*X(mind))
hX = sln_AISoftImpute1 * maxS;

meaMat(i,:) = AE(hX, X, Xmis, group);
% NAE(sln_AISoftImpute1, X / maxS, Xmis/maxS, group)
end

mean(meaMat)
save ./simu_R1/simuData_Mat/Mdata_R1/AE_simu3_misRate2mechMAR_PLAIS qVec meaMat
