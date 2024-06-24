function [ success_or_not ] = SimpleShootingStiefelChecks( Delta_rec, Y0, Y1, Delta_exact, tol )

% function [ success_or_not ] = SimpleShootingStiefelChecks( Delta_k, Y0, Y1, Delta_exact, tol )
% Purpose: Contains all the checks to be performed on the solution Delta_k
%          found by the Simple Shooting on the Stiefel manifold.
% Created:     19.09.2016
% Last change: 09.09.2020

%   Mar 22, 2024:
%   September 9, 2020:
%       Added 'Norm_Delta_exact' and 'Norm_Delta_reconstructed'.

[ n, p ] = size(Delta_rec);

disp('+-------------------------------------------------+')
disp('|             CHECKS SINGLE SHOOTING              |')
disp('+-------------------------------------------------+')

% % 1) Check that Fs is identically zero (component-wise)
% if abs( FDelta(end) ) < tol
%     disp('F(Delta) is zero:                   OK')
% else
%     disp('F(Delta) is zero:                   NO :-(')
% end

% % 2) Is the Delta_k the same as Delta_exact?
% if abs( Delta_k - Delta_exact ) < tol*ones(n,p)
%     disp('Delta_k is equal to Delta_exact:    OK')
% else
%     disp('Delta_k is equal to Delta_exact:    NO')
% end

% 3) Is Y1 reached by Delta_k?
if abs( Stiefel_Exp( Y0, Delta_rec ) - Y1 ) < tol*ones(n,p)
    disp('Y1 is reached by Delta_k:                       OK')
    success_or_not = 1;
else
    disp('Y1 is reached by Delta_k:                       NO :-(')
    success_or_not = 0;
end

Norm_Delta_exact = GetCanonicalNormDelta( Delta_exact, Y0 );
Norm_Delta_rec = GetCanonicalNormDelta( Delta_rec, Y0 );

% 4) Is the Norm_Delta_exact the same as Norm_Delta_reconstructed?
if abs( Norm_Delta_exact - Norm_Delta_rec ) < tol
    disp('Norm_Delta_exact is equal to Norm_Delta_rec:    OK')
else
    disp('Norm_Delta_exact is equal to Norm_Delta_rec:    NO')
end
%     formatSpec = ' %2d   %10.4e   %10.4e   %10.4f  \n';
fprintf( '   Norm_Delta_exact: %10.4f \n', Norm_Delta_exact )
fprintf( '   Norm_Delta_rec:   %10.4f \n', Norm_Delta_rec )

end