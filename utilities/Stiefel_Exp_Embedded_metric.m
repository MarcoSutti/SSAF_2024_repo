function [ Y ] = Stiefel_Exp_Embedded_metric( X, Delta )

% function [ Y ] = Stiefel_Exp_Embedded_metric( X, Delta )
% Purpose: Compute the Riemannian exponential map on the Stiefel manifold
%          with the embedded Euclidean metric.
% Reference: Eq. (5.26) in Absil, Mahony, and Sepulchre. Or eq. (6) in
%            Bryner's paper. Or eq. (6) in Zimmermann and Huper.
% Created:     2024.01.21
% Last change: 2024.02.20

%   Jan 21, 2024:
%       Created.

[ ~, p ] = size(X);

A0 = X'*Delta;
S0 = Delta'*Delta;

expm_mat = expm( [ A0, -S0; eye(p), A0 ] );

% Exponential mapping
Y = [ X, Delta ] * (expm_mat(:,1:p) * expm(-A0));

% Eq. (7) in Zimmermann and Huper:
% Y = expm(Delta*X'-X*Delta')*X*expm(-A0);

end