function [ Y0_tilde, Y1_tilde, Y0perp_tilde, U1 ] = SetBabyProblem( Y0, Y1 )

% function [ Y0_tilde, Y1_tilde, Y0perp_tilde, U1 ] = SetBabyProblem( n, p, Y0, Y1 )
% Created:     14.11.2016
% Last change: 14.11.2016

p = size( Y0, 2 );

M = Y0'*Y1;

[ U1, N ] = qrPosDiagR( Y1 - Y0*M );  % 2023.08.08: This is also a kind of bottleneck...

% We solve a geodesic problem on St(2*p,p) with the following initial
% conditions:
Y1_tilde = [ M; N ];
Ipp = speye(2*p);
Y0_tilde = Ipp( :, 1:p );
Y0perp_tilde = Ipp( : , p+1:end );

end

