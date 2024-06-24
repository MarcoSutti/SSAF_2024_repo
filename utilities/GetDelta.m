function [ Delta_n ] = GetDelta( Y0 )

% function [ Delta_n ] = GetDelta( Y0 )
% Purpose: Returns a random, unit-norm tangent vector in the tangent space
%          T_{Y0}St(n,p).
% Last change: 2016.07.18
% Last change: 2024.01.18

[ n, p ] = size(Y0);

Omega   = rand( p, p );
Omega   = Omega - Omega';          % random, p-by-p skew-symmetric matrix
K_0     = rand( n, p );         % random, free n-by-p matrix
Delta = Y0*Omega + K_0-Y0*(Y0'*K_0);

% Normalize Delta w.r.t. the "canonical" metric on the tangent space
norm_Delta = GetCanonicalNormDelta( Delta, Y0 );
Delta_n = Delta/norm_Delta;

end