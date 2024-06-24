%==========================================================================
% Computes the matrix exponential expm of a unit norm skew-symmetric
% matrix in R1000×1000 is computed in 0.45 seconds (time is averaged over
% 100 runs).
% We give this reference for the computational cost of the matrix
% exponential expm following [13, Section 5.2].

% Created:     2024.02.23
% Last change: 2024.06.24
%==========================================================================

close all; clear; clc;

A = rand(1000); Omega = .5*(A-A'); Omega = Omega/norm(Omega);


nb_experiments = 100;

for idx_experiment=1:nb_experiments
    tic
    expm(A);
    time(idx_experiment) = toc;
end


fprintf( '  Number of random experiments:   %d. \n', nb_experiments );
fprintf( '  Averaged computational time :   %5.5f. \n', mean(time) );
