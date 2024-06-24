%==========================================================================
% Driver for the Simple Shooting Method on the Stiefel Manifold.
% 2024: New algorithm with approximation of the Fréchet derivative!!!!

% Created:     2024.02.17
% Last change: 2024.06.24

%   Feb 17, 2024:
%       Created.
%==========================================================================

% Sets default graphics interpreter, paths and colors.
SSLF_startup;

% rng(33)

%--------------------------------------------------------------------------
% Data
%--------------------------------------------------------------------------
% Set dimensions of St(n,p)
n = 1500;
p = 500;

% Fix a distance distY0Y1 between the endpoints Y0 and Y1
prefactor_distY0Y1 = 0.5;
distY0Y1 = prefactor_distY0Y1 * pi;   % Scegliamo un multiplo di pi greco.

param.tol = 1e-5;
param.maxiter = 100;
param.verbose = 0;

nb_experiments = 10;

time_SAF_array = zeros(1,nb_experiments);
nb_its_SAF  = zeros(1,nb_experiments);
time_SEJ_array = zeros(1,nb_experiments);
nb_its_SEJ  = zeros(1,nb_experiments);

%--------------------------------------------------------------------------
% END OF DATA
%--------------------------------------------------------------------------

for idx_experiment = 1:nb_experiments

    % Create random Stiefel matrix Y0
    Y0 = orth( rand( n, p ) );
    Y0perp = null(Y0');    % The columns of Y0perp span the orthogonal complement to the subspace span(Y0)
    Delta_exact = distY0Y1 * GetDelta( Y0 );

    % Map the tg vector onto the manifold
    [ Y1 ] = Stiefel_Exp( Y0, Delta_exact );

    % rank( Y1 - Y0*(Y0'*Y1), 1e-10 )

    [ iter_SAF, time_SAF, Delta_rec_SAF, ~ ] = SingleShooting_ApproxFrechet( Y0, Y1, Y0perp, param );
    time_SAF_array(idx_experiment) = time_SAF;
    nb_its_SAF(idx_experiment) = iter_SAF;
    % [ success_SAF ] = SimpleShootingStiefelChecks( Delta_rec_SAF, Y0, Y1, Delta_exact, param.tol );

end

if nb_experiments > 1
    fprintf( '  Number of random experiments:   %d. \n', nb_experiments );
    disp('+-------------------------------------------------+')
    disp('| SHOOTING WITH APPROX OF THE FRECHET DERIVATIVE  |')
    disp('+-------------------------------------------------+')
    fprintf( '  Averaged computational time :   %5.5f. \n', mean(time_SAF_array) );
    fprintf( '  Averaged nb of iterations   :   %5.2f. \n', mean(nb_its_SAF) );
    % fprintf( '  Success   :   %5.2f. \n', mean(success_SAF)*100 );
end

%--------------------------------------------------------------------------
% Plot.
%--------------------------------------------------------------------------

% PlotConvergenceSingleShooting( norm_update )
