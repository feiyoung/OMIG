function [measure] = measurefun(X, hX, type)
% function to evaluate the smallest nonzero canonical correlation between
% two set of variables.
%----------------------------------------------------------------------------------
% Author:       Liu Wei
% Maintainer:    Liu Wei <weiliu@smail.swufe.edu.cn>
% Date:         May. 25, 2020
% Copyright (c) 2020, Liu Wei
% All rights reserved.

if(~exist('type', 'var'))
    type = 'canonical';
end
q = size(X, 2);
switch type 
    case 'canonical'
        [~, ~, r] = canoncorr(hX,X);
        measure = r(end);
    case 'Fnorm'
        measure = norm(X- hX, 'fro')^2/ prod(size(X));
end