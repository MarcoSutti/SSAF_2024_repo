function [ iter, norm_update, Delta_k, param ] = SimpleShootingStiefel_ApproxFrechet( Y0, Y1, Y0perp, Delta_0, param )

% function [ iter, norm_update, Delta_k, param ] = SimpleShootingStiefel_ApproxFrechet( Y0, Y1, Y0perp, Delta_0, param )
% Encapsulates all that is needed to do the simple shooting with the
% formulation using the parametrization of tangent vectors.
% Created:     2024.02.17
% Last change: 2024.03.15

if param.verbose==2
    fprintf( ' it    residual\n' );
end

[ n, p ] = size(Y0);

% Muovi questi due calcoli fuori, sono sempre uguali!!!
Q = [ Y0, Y0perp ];

% Q'*Y1 it is done once and for all here:
QtY1 = Q'*Y1;

% Initializations:
iter = 1;
param.flag  = true;   % we initialize the flag to true
norm_update(1) = param.tol + 1;
Delta_k = Delta_0;
Onminp =  zeros( n-p );

while and( norm_update(iter) > param.tol, iter < param.maxiter )

    iter = iter + 1;

    Omega = Y0'*Delta_k;
    K = Y0perp'*Delta_k;

    hQ = [ Omega,   -K';
               K,    Onminp ];

    % check = full(Y0perp*Y0perp')

    exphQ = expm(hQ);

    %     norm(.5*Y0tDelta_k,"fro")

    % MS, 2024.02.26:
    % Approximation of the Fréchet derivative of the matrix exponential
    RHS = QtY1 - exphQ(:,1:p);

    KtK = K'*K;
    A = speye(p) + .5*Omega + .25*KtK;
    B = .5*Omega - .25*KtK;
    N = RHS(p+1:n,:);
    
    M = speye(p) + .5*Omega;

    KtN = K'*N;
    KtNinvM = KtN/M;

    C = RHS(1:p,:) + .5*(KtNinvM+KtNinvM');

    dOmega = lyap(A,B,-C);

    dK = (N - .5*K*dOmega)/M;

%     dOmega = RHS(1:p,:);
%     dK = RHS(p+1:n,:);

    dDelta_k = Y0*dOmega + Y0perp * dK;

    % NO line-search at the moment.
%     fun = @(t) eval_F_LS( QtY1, Onminp, Y0, Y0perp, p, Delta_k, dDelta_k, t );
% 
%     t_star = fminbnd( fun, 1e-8, 10 );
% 
% %     if param.verbose==2
%         formatSpec = '      SS-LS, t_star:  %1.4f\n';
%         fprintf( formatSpec, t_star )
% %     end

%     if iter < 100
%         t_star = 0.5;
%     else
%         t_star = 1;
%     end

    % Use unitary step size
%     t_star = 1;

    Delta_k = Delta_k + dDelta_k;

%     Delta_k = Y0*(Omega + dOmega) + Y0perp * (K + dK);
    Delta_k = ProjTgSpaceStiefel( Y0, Delta_k );

    norm_update(iter) = norm( dDelta_k, "fro" );

    if param.verbose==2
        formatSpec = ' %2d   %10.4e\n';
        fprintf( formatSpec, iter, norm_update(iter) )
    end

end

if or( iter==param.maxiter, or( norm_update(iter) > param.tol, isnan(norm_update(iter)) ) )
    % 18.12.2016: if we reach the maximum number of iterations allowed for
    % Single Shooting, we set the flag to 'false'
    if iter==param.maxiter
        disp( 'SINGLE SHOOTING: REACHED MAXIMUM NUMBER OF ITERATIONS!!!' )
    end
    if or( norm_update(iter) > param.tol, isnan(norm_update(iter)) )
        disp( 'SINGLE SHOOTING: norm_update BLOWED UP!!!' )
    end
    param.flag = false;
    Delta_k = NaN;
    return;
end

end


% function [ res ] = eval_F_LS( QtY1, Onminp, Y0, Y0perp, p, Delta_k, dDelta_k, alpha )
% 
% Delta_k = Delta_k + alpha * dDelta_k;
% Delta_k = ProjTgSpaceStiefel( Y0, Delta_k );
% 
% Omega_k = Y0'*Delta_k;
% K_k = Y0perp'*Delta_k;
% 
% A = [ Omega_k,     -K_k';
%     K_k,    Onminp ];
% 
% expmA = expm(A);
% 
% res = norm( QtY1 - expmA(:,1:p), "fro" );
% 
% end

