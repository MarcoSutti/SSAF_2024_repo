function [ PXZ ] = ProjTgSpaceStiefel( X, Z )

% function [ PXZ ] = ProjTgSpaceStiefel( X, Z )
% Purpose: Projects Z onto the tangent space T_{X}St(n,p).
% Created:     08.08.2016
% Last change: 26.08.2016

% Z = reshape(Z,[n,p]);
% X = reshape(X,[n,p]);
XtZ = X'*Z;

PXZ = Z - (1/2)*X*( XtZ + XtZ' );
%PXZ = PXZ(:);

end