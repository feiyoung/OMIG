function [hH] = localupdateH2(Xmis, gcell, hBm, type, dropout)
% function to update latent factor matrix H.
%----------------------------------------------------------------------------------
% Author:       Liu Wei
% Maintainer:    Liu Wei <weiliu@smail.swufe.edu.cn>
% Date:         May. 25, 2019
% Copyright (c) 2019, Liu Wei
% All rights reserved.

mind = ~isnan(Xmis);
n = size(Xmis, 1);
q_1 = size(hBm, 2);
ng = size(type,1);
if all(dropout ~=0) && (~isempty(setdiff(dropout , 1:ng)))
    error('dropout setting is wrong!')
end
idres = setdiff(1:ng, dropout);
Harray = zeros(n,q_1, length(idres));
w = ones(1,n);
for j = idres
    H2 = zeros(q_1,n);
    if strcmp(type{j,1}, 'normal')
        
        parfor i = 1:n
            ij_obs = intersect(gcell{j}, find(mind(i,:)==1));
            w = 1./ (std(Xmis(:, ij_obs), 'omitnan').^2);
            H2(:,i) = glmfit(hBm(ij_obs,:), Xmis(i,ij_obs)',type{j,1},'link',type{j,2}, 'constant', 'off', 'weights', w);
        
        end
    elseif strcmp(type{j,1}, 'binomial')
        p1 = length(gcell{j});
        parfor i = 1:n
            ij_obs = intersect(gcell{j}, find(mind(i,:)==1));
            ntrail_i = length(unique(Xmis(i,ij_obs)))-1;
           
            H2(:,i) = glmfit(hBm(ij_obs,:), [Xmis(i,ij_obs)',ntrail_i*ones(length(ij_obs),1)],type{j,1},'link',type{j,2}, 'constant', 'off');
            
        end
    else
        parfor i = 1:n
            ij_obs = intersect(gcell{j}, find(mind(i,:)==1));
            H2(:,i) = glmfit(hBm(gcell{j},:), Xmis(i,gcell{j})',type{j,1},'link',type{j,2}, 'constant', 'off');
            
        end
    end
    
    Harray(:,:,j) = H2';
end
hH = mean(Harray, 3);