function [ iter, norm_update, Delta_rec, param ] = SimpleShootingStiefel_Baby_ApproxFrechet( Y0, Y1, param )

% function [ iter, norm_update, Delta_rec, param ] = SimpleShootingStiefel_Baby_ApproxFrechet( Y0, Y1, param )
% Created:     2024.02.17
% Last change: 2024.02.17

[ Y0_tilde, Y1_tilde, Y0perp_tilde, U1 ] = SetBabyProblem( Y0, Y1 );

Delta_0 = GetStartingGuessDelta( Y0_tilde, Y1_tilde );

[ iter, norm_update, Delta_k, param ] = SimpleShootingStiefel_ApproxFrechet( Y0_tilde, Y1_tilde, Y0perp_tilde, Delta_0, param );

% We reconstruct the tangent vector from this solution
A_tilde = Y0_tilde'*Delta_k;
B_tilde = Y0perp_tilde'*Delta_k;
Delta_rec = Y0*A_tilde + U1*B_tilde;

end