function [ k, Delta, err_Bryner ] = Bryner_Stiefel( Yinit, Yfin, param_Bryner )

% function [ k, Delta, time_Bryner ] = Bryner_Stiefel( Yinit, Yfin, embedded_norm_Delta_exact, param_Bryner )
% Purpose: runs Bryner's method on the Stiefel manifold.
% Reference: Zimmermann and Huper, 2022, Algorithm 1.

% Created:     2024.01.21
% Last change: 2024.02.23

%   Jan 24, 2024:
%       Removed plot of wrong error quantity
%   Jan 21, 2024:
%       Created.

[ n, p ] = size(Yinit);

if param_Bryner.verbose>1
    disp('+-------------------------------------------------+')
    disp('|              BRYNERS SHOOTING METHOD            |')
    formatSpec = '|                     T = %1.0d                      |\n';
    fprintf( formatSpec, param_Bryner.T );
    disp('+-------------------------------------------------+')
end

k = 1;
t = linspace(0,1,param_Bryner.T);
gamma = norm( Yinit - Yfin, "fro");

err_Bryner(1) = gamma;

Proj = ProjTgSpaceStiefel( Yinit, Yfin-Yinit );
Delta = (gamma/norm(Proj,"fro"))*Proj;
% norm_error_L_omega_Bryner(1) = abs( norm( Delta, "fro" ) - embedded_norm_Delta_exact );

Uvec = zeros( n*p, param_Bryner.T );

while (err_Bryner(k)>param_Bryner.tol && k<param_Bryner.maxit)
    if param_Bryner.verbose>1
        formatSpec = ' %3d  |U(1)-U_1| %0.4e\n';
        fprintf( formatSpec, k, err_Bryner(k))
    end

    for j=1:param_Bryner.T
        % Shot discrete geodesic with tangent vector Delta.
        Uvec(:,j) = reshape( Stiefel_Exp_Embedded_metric( Yinit, t(j)*Delta ), [ n*p, 1 ] );
    end
    Delta_s = reshape( Uvec(:,end) - Yfin(:), [n, p] );
	gamma = norm( Delta_s, "fro" );
	err_Bryner(k+1) = gamma;

	for j=param_Bryner.T:-1:1        % kind of backpropagation using parallel transport.
		Proj = ProjTgSpaceStiefel( reshape( Uvec(:,j), [n, p] ), Delta_s );
		Delta_s = (gamma/norm(Proj,"fro"))*Proj;
	end
	
    Delta = Delta - Delta_s;   % no line search? Vedi ultima frase pagina 1145.
%     norm_error_L_omega_Bryner(k+1) = abs( norm( Delta, "fro" ) - embedded_norm_Delta_exact );%% Tom: ATTENTO QUI! Un vettore tangente pressupone di usare Exp_embedded (delta), l'atro invece no, e' solamente normalizzato w.r.t. un'altra metrica!
    k = k+1;
end

if param_Bryner.verbose>1
    %formatSpec = ' %3d    %0.4e\n';
    formatSpec = ' %3d  |U(1)-U_1| %0.4e\n';
    fprintf( formatSpec, k, err_Bryner(k))
end

if param_Bryner.verbose>1
    formatSpec = '  Bryners shooting converged in %d iterations.\n';
    fprintf( formatSpec, k );
end

% if param_Bryner.verbose>0 %DA CANCELLARE
%     fprintf( "   Total time Bryner's method: %2.5f s.\n", time_Bryner );
% end
end