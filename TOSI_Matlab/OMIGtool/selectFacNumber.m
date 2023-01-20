function [icMat] = selectFacNumber(Xmis,group, type, q_set,lambda,parallel, verbose)
% function to evluate the information criteria values for selecting factor
% number.
%----------------------------------------------------------------------------------
% Author:       Liu Wei
% Maintainer:    Liu Wei <weiliu@smail.swufe.edu.cn>
% Date:         June. 25, 2021
% Copyright (c) 2021, Liu Wei
% All rights reserved.
if(~exist('verbose', 'var') || isempty(verbose))
    verbose = 0;
end
nq_set = length(q_set);
icMat = zeros(nq_set, 2);

 for i = 1: nq_set
        [~,~,~,~,icMat(i,:)] = OrMIG(Xmis, group, type, q_set(i), 1e-5, 10, lambda, verbose, parallel);
 end

end
