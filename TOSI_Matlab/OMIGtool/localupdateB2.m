function [B1] = localupdateB2(Xmis1, hH, type1,  parallel)
% function to update loading matrix B.
%----------------------------------------------------------------------------------
% Author:       Liu Wei
% Maintainer:    Liu Wei <weiliu@smail.swufe.edu.cn>
% Date:         May. 25, 2021
% Copyright (c) 2021, Liu Wei
% All rights reserved.
% Xmis1 = Xmis(:,gcell{jj});
% type1 = type(jj,:)
n = size(Xmis1, 1);
q = size(hH,2);
p1 = size(Xmis1, 2);
B1 = zeros(q+1, p1);
jg = 1:p1;

if strcmp(type1{1}, 'binomial')
    if parallel
        parfor j = jg
            
            B1(:,j) = glmfit_binomial(hH, Xmis1(:,j), type1);
        end
    else
        for j = jg
            
            B1(:,j) = glmfit_binomial(hH, Xmis1(:,j), type1);
        end
    end
else
    if parallel
        parfor j = jg
            B1(:,j) = glmfit_othertype(hH, Xmis1(:,j), type1);
        end
    else
        for j = jg
            B1(:,j) = glmfit_othertype(hH, Xmis1(:,j), type1);
        end
    end
end


function [beta] = glmfit_othertype(hH, Xj, type1)
ind = ~isnan(Xj);
beta = glmfit(hH(ind==1,:), Xj(ind==1),type1{1},'link',type1{2}, 'constant', 'on');

function [beta] = glmfit_binomial(hH, Xj, type1)
ind = ~isnan(Xj);
ntrail_j = length(unique(Xj(ind==1)))-1;
beta = glmfit(hH(ind==1,:), [Xj(ind==1,:), ntrail_j*ones(sum(ind),1)],type1{1},'link',type1{2}, 'constant', 'on');


        
