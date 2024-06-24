function [ ] = BrynerShootingChecks( Delta_Bryner, Y0, Y1, param_Bryner )

% function [ ] = BrynerShootingChecks( Delta_Bryner, Y0, Y1, Delta_exact, param_Bryner )
% Purpose: Contains all the checks to be performed on the solution Delta_k
%          found by the Simple Shooting on the Stiefel manifold.
% Created:     2024.02.21
% Last change: 2024.02.21

%   Feb 21, 2020:
%       Created.

[ n, p ] = size(Delta_Bryner);

disp('+-------------------------------------------------+')
disp('|            CHECKS BRYNERS SHOOTING              |')
disp('+-------------------------------------------------+')

% Check norm of the difference Delta_exact and Delta_rec_Bryner:
% check_norm = norm( Delta_exact - Delta_Bryner, "fro" );
% fprintf( "   Norm(Delta_exact - Delta_reconstructed_Bryner) = %.4e.\n", check_norm );

% % 2) Is the Delta_k the same as Delta_exact?
% if abs( Delta_exact - Delta_Bryner ) < param_Bryner.tol*ones(n,p)
%     disp('Delta_k is equal to Delta_exact:    OK')
% else
%     disp('Delta_k is equal to Delta_exact:    NO')
% end

% 3) Is Yfin reached by Delta_Bryner?
check_norm = norm( Stiefel_Exp_Embedded_metric( Y0, Delta_Bryner ) - Y1, "fro" );
fprintf( "   norm( Stiefel_Exp_Embedded_metric( Y0, Delta_Bryner ) - Y1 ) = %.4e.\n", check_norm );

if abs( Stiefel_Exp_Embedded_metric( Y0, Delta_Bryner ) - Y1 ) < param_Bryner.tol*ones(n,p)
    disp('Yfin is reached by Delta_Bryner:    OK')
else
    disp('Yfin is reached by Delta_Bryner:    NO :-(')
end

% Norm_Delta_exact = norm( Delta_exact, "fro" );
Norm_Delta_Bryner = norm( Delta_Bryner, "fro" );

%     formatSpec = ' %2d   %10.4e   %10.4e   %10.4f  \n';
% fprintf( '   Norm_Delta_exact:         %10.4f \n', Norm_Delta_exact )
fprintf( '   Norm_Delta_Bryner:        %10.4f \n', Norm_Delta_Bryner )

end