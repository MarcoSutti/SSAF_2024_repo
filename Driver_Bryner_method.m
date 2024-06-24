%==========================================================================
% Driver for Bryner's method for computing geodesics on the Stiefel
% manifold.

% Created:     2024.02.20
% Last change: 2024.02.20

%   Feb 20, 2024:
%       Created.
%==========================================================================

close all; clear; clc;

addpath(genpath('.'))

options_plot;

%--------------------------------------------------------------------------
param_Bryner.T = 20;   % This is the parameter T in the paper. It is a
% kind of resolution for approximating the parallel translation.

% Stopping criterion parameters for Bryner's method:
param_Bryner.maxit = 500;   % max number of iterations
param_Bryner.tol = 1e-3;    % tolerance

param_Bryner.verbose = 2;
%--------------------------------------------------------------------------


%--------------------------------------------------------------------------
% Data
%--------------------------------------------------------------------------
% Set dimensions of St(n,p)
n = 12;
p = 3;

% Fix stream of random numbers for reproducibility
% rng(23)

% Fix a distance between the endpoints X and Y:
prefactor_distXY = 1.30;
distXY = prefactor_distXY * pi;   % Scegliamo un multiplo di pi greco.

%--------------------------------------------------------------------------
% End of data
%--------------------------------------------------------------------------

% Create random Stiefel matrix X
Yinit = orth( rand( n, p ) );

% Create a random tangent vector Delta in T_{X}St(n,p)
Delta_exact = distXY * GetDelta( Yinit );

% Map the tg vector onto the manifold:
[ Yfin ] = Stiefel_Exp_Embedded_metric( Yinit, Delta_exact );

% canonical_norm_Delta_exact = GetCanonicalNormDelta( Delta_exact, Yinit );
% embedded_norm_Delta_exact = norm( Delta_exact, "fro" );

%--------------------------------------------------------------------------
% BRYNER'S METHOD ON THE STIEFEL MANIFOLD
%--------------------------------------------------------------------------
% Run Bryner's algorithm:
tic;
[ nb_its, Delta_Bryner, err_L_Bryner ] = Bryner_Stiefel( Yinit, Yfin, param_Bryner );
time_Bryner = toc;
formatSpec = '  Computational time:   %5.5f. \n';
fprintf( formatSpec, time_Bryner );
%--------------------------------------------------------------------------

%--------------------------------------------------------------------------
% Postprocessing of Bryner's method
%--------------------------------------------------------------------------
BrynerShootingChecks( Delta_Bryner, Yinit, Yfin, param_Bryner );


%--------------------------------------------------------------------------
% Convergence plot.
%--------------------------------------------------------------------------
semilogy( 1:length(err_L_Bryner), err_L_Bryner, ...
    ':', 'Color', blue, 'LineWidth', 3 );
hold on;
grid on;
xlabel('Number of iterations', 'FontSize', 18 )
ylabel('Error', 'FontSize', 18 )

% Legend
xlabel('Number of iterations', 'FontSize', 14 )
ylabel('Error', 'FontSize', 14 )
