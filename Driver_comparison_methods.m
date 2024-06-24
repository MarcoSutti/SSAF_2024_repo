%==========================================================================
% Driver for comparing my Single Shooting method with Bryner's shooting
% method.

% Created:     2024.02.20
% Last change: 2024.04.05

%   Mar 20, 2024:
%       Added Zimmermann's 2017 method from:
%          SIAM J. Matrix Anal. Appl. Vol. 38, No. 2, pp. 322–342.
%   Feb 20, 2024:
%       Created.
%==========================================================================

% Sets default graphics interpreter, paths and colors.
SSLF_startup;

options_plot;

%--------------------------------------------------------------------------
% Data
%--------------------------------------------------------------------------
% Set dimensions of St(n,p)
n = 500;   % 2000, 120
p = 128;    % 500, 30

% Fix a length distXY
prefactor_distXY = 0.5;
distXY = prefactor_distXY * pi;   % Scegliamo un multiplo di pi greco.

%--------------------------------------------------------------------------
param_Bryner.T = 20;   % This is the parameter T in the paper. It is a
% kind of resolution for approximating the parallel translation.

% Stopping criterion parameters for Bryner's method:
param_Bryner.maxit = 100;   % max number of iterations
param_Bryner.tol = 1e-3;    % tolerance

param_Bryner.verbose = 1;
%--------------------------------------------------------------------------
% Parameters Single Shooting
param_SS.tol = 1e-3;
param_SS.maxiter = 500;

param_SS.verbose = 1;
%--------------------------------------------------------------------------
% Parameters for Leap-Frog Single Shooting
param_LF.tolLF = 1e-3;
param_LF.maxiterLF = 100;
param_LF.tol = 1e-3;
param_LF.maxiter = 100;

param_LF.verboseLF = 0;
param_LF.verbose = 0;

% Number of subintervals for leapfrog
param_LF.m = 5;
%--------------------------------------------------------------------------
% Number of repeated random experiments
nb_experiments = 100;
%--------------------------------------------------------------------------
% END OF DATA
%--------------------------------------------------------------------------

time_Bryner = zeros(1,nb_experiments);
nb_its_Bryner  = zeros(1,nb_experiments);
time_SS_array = zeros(1,nb_experiments);
nb_its_SS  = zeros(1,nb_experiments);
time_Zimmermann = zeros(1,nb_experiments);
nb_its_Zimmermann  = zeros(1,nb_experiments);
time_LF_array = zeros(1,nb_experiments);
nb_its_LF  = zeros(1,nb_experiments);

for idx_experiment=1:nb_experiments

    % Create random Stiefel matrix Y0
    Y0 = orth( rand( n, p ) );
    Y0perp = null(Y0');    % The columns of Y0perp span the orthogonal complement to the subspace span(Y0)
    % Create a random tangent vector Delta in T_{Y0}St(n,p) with a given length
    % distXY
    Delta_exact = distXY * GetDelta( Y0 );

    % Map the tg vector onto the manifold
    [ Y1 ] = Stiefel_Exp( Y0, Delta_exact );

    %--------------------------------------------------------------------------
    % BRYNER'S METHOD ON THE STIEFEL MANIFOLD
    %--------------------------------------------------------------------------
    % Run Bryner's algorithm:
    tic;
    [ iter_Bryner, Delta_Bryner, err_L_Bryner ] = Bryner_Stiefel( Y0, Y1, param_Bryner );
    time_Bryner(idx_experiment) = toc;
    nb_its_Bryner(idx_experiment) = iter_Bryner;

    %--------------------------------------------------------------------------
    % ZIMMERMANN'S 2017 METHOD
    %--------------------------------------------------------------------------
    tic
    [ Delta_Zimmermann, iter_Zimmermann, conv_hist, norm_logV0 ] = Stiefel_Log_supp(Y0, Y1, param_SS.tol);
    time_Zimmermann(idx_experiment) = toc;
    nb_its_Zimmermann(idx_experiment) = iter_Zimmermann;
%     disp('+-------------------------------------------------+')
%     disp('|              ZIMMERMANN 2017 METHOD             |')
%     disp('+-------------------------------------------------+')
%     success_Zimmermann(idx_experiment) = SimpleShootingStiefelChecks( Delta_Zimmermann, Y0, Y1, Delta_exact, param_SS.tol );

    %--------------------------------------------------------------------------
    % SHOOTING METHOD WITH APPROXIMATION OF THE JACOBIAN
    %--------------------------------------------------------------------------
    [ iter_SS, time_SS, Delta_SS, norm_update ] = SingleShooting_ApproxFrechet( Y0, Y1, Y0perp, param_SS );
    time_SS_array(idx_experiment) = time_SS;
    nb_its_SS(idx_experiment) = iter_SS;
%     disp('+-------------------------------------------------+')
%     disp('| SHOOTING WITH APPROX OF THE FRECHET DERIVATIVE  |')
%     disp('+-------------------------------------------------+')
%     success_SS(idx_experiment) = SimpleShootingStiefelChecks( Delta_SS, Y0, Y1, Delta_exact, param_SS.tol );

end

if nb_experiments > 1
    fprintf( '  Number of random experiments:   %d. \n', nb_experiments );
    disp('+-------------------------------------------------+')
    disp('|              BRYNERS SHOOTING METHOD            |')
    formatSpec = '|                     T = %1.0d                      |\n';
    fprintf( formatSpec, param_Bryner.T );
    disp('+-------------------------------------------------+')
    fprintf( '  Averaged computational time :   %5.5f. \n', mean(time_Bryner) );
    fprintf( '  Averaged nb of iterations   :   %5.2f. \n', mean(nb_its_Bryner) );
    disp('+-------------------------------------------------+')
    disp('|              ZIMMERMANN 2017 METHOD             |')
    disp('+-------------------------------------------------+')
    fprintf( '  Averaged computational time :   %5.5f. \n', mean(time_Zimmermann) );
    fprintf( '  Averaged nb of iterations   :   %5.2f. \n', mean(nb_its_Zimmermann) );
%     fprintf( '  Success   :   %5.2f. \n', mean(success_Zimmermann)*100 );
    disp('+-------------------------------------------------+')
    disp('|   SINGLE SHOOTING WITH APPROX OF THE JACOBIAN   |')
    disp('+-------------------------------------------------+')
    fprintf( '  Averaged computational time :   %5.5f. \n', mean(time_SS_array) );
    fprintf( '  Averaged nb of iterations   :   %5.2f. \n', mean(nb_its_SS) );
%     fprintf( '  Success   :   %5.2f. \n', mean(success_SS)*100 );
end


%--------------------------------------------------------------------------
% Postprocessing of Bryner's method
%--------------------------------------------------------------------------
% BrynerShootingChecks( Delta_Bryner, Y0, Y1, param_Bryner );


%--------------------------------------------------------------------------
% Postprocessing of Single Shooting
%--------------------------------------------------------------------------
% All the checks:
% SimpleShootingStiefelChecks( Delta_SS, Y0, Y1, Delta_exact, param_SS.tol )
