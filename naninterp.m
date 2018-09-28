function X = naninterp(X)
% This function uses cubic interpolation to replace NaNs
% See INTERP1 for more info
X(isnan(X)) = interp1(find(~isnan(X)), X(~isnan(X)), find(isnan(X)),'cubic');
