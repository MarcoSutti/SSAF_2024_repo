function [ iter, time_SS, Delta_rec, norm_update ] = SingleShooting_ApproxFrechet( Y0, Y1, Y0perp, param )

% function [ Delta_rec, norm_update ] = SingleShooting_ApproxFrechet( Y0, Y1, param )
% Purpose:
% Created:     2024.02.21
% Last change: 2024.02.23

%   Feb 21, 2024:
%       Created.

[ n, p ] = size(Y0);

if param.verbose>1
    disp('+-------------------------------------------------+')
    disp('| SHOOTING WITH APPROX OF THE FRECHET DERIVATIVE  |')
    disp('+-------------------------------------------------+')
end

if p < n/2
    tic;
    [ iter, norm_update, Delta_rec, param ] = SimpleShootingStiefel_Baby_ApproxFrechet( Y0, Y1, param );
    time_SS = toc;
else
    tic;
    Delta_0 = GetStartingGuessDelta( Y0, Y1 );
    [ iter, norm_update, Delta_rec, param ] = SimpleShootingStiefel_ApproxFrechet( Y0, Y1, Y0perp, Delta_0, param );
    time_SS = toc;
end

if param.verbose>1
    formatSpec = '  Computational time:   %5.5f. \n';
    fprintf( formatSpec, time_SS );
    if param.flag==1
        formatSpec = '  Single shooting converged in %d iterations.\n';
        fprintf( formatSpec, iter );
    else
        disp( 'Single shooting failed.' )
    end
end

end