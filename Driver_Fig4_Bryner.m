%==========================================================================
% Driver for plotting comparison between my Single Shooting method
% with Bryner's shooting method.

% Created:     2024.02.23
% Last change: 2024.06.24

%   Apr 5, 2024:
%       Added new values of avg. computational times.
%   Feb 23, 2024:
%       Created.
%==========================================================================

% Sets default graphics interpreter, paths and colors.
SSLF_startup;

options_plot;

color_array = [ blue; red; green; yellow; Cerulean ];
style_array = {':'; '--'; '-'; '-.'; ':' };

linewidth = 2.5;
marker_size = 10;

% LineWidth of the MarkerEdge:
myMarkerLineWidth = 0.5;
%--------------------------------------------------------------------------

% Doubling values of n, from 10 to 10240:
n_array = 10*2.^(0:10);

% Bryners_times = [ 0.00333, 0.00289, 0.00259, 0.00260, 0.00248, 0.00257, ...
%     0.00288, 0.00412, 0.00741, 0.01066, 0.01751 ];
% 
% my_shooting_times = [ 0.00184, 0.00166, 0.00149, 0.00149, 0.00144, ...
%     0.00149, 0.00231, 0.00248, 0.00232  ];

Bryners_times_vs_n = [ 0.00400, 0.00367, 0.00337, 0.00312, 0.00310, 0.00328, ...
    0.00371, 0.00543, 0.00856, 0.01056, 0.01596 ];

Zimmermanns_times_vs_n = [ 0.00091, 0.00093, 0.00095, 0.00101, 0.00105, 0.00107, ...
    0.00105, 0.00104, 0.00135, 0.00131, 0.00144 ];

SSAF_times_vs_n = [ 0.00080, 0.00091, 0.00075, 0.00081, 0.00086, 0.00096, ...
    0.00091, 0.00100, 0.00121, 0.00132, 0.00141];

%--------------------------------------------------------------------------

%--------------------------------------------------------------------------
% PLOT
%--------------------------------------------------------------------------
figure( 'units', 'normalized', 'outerposition', [ 0 1 .66 .5 ] );

tiledlayout( 1, 2 );
% tiledlayout( 1, 3 , 'Padding', 'tight' );
% tiledlayout( 1, 3, 'TileSpacing', 'compact', 'Padding', 'compact');

nexttile
handle_plot_vs_n(1) = loglog( n_array, Bryners_times_vs_n, 'o', ...
    'LineStyle', style_array(1), 'Color', color_array(1,:), 'LineWidth', linewidth, ...
    'MarkerEdgeColor', color_array(1,:), 'MarkerFaceColor', color_array(1,:), ...
    'MarkerSize', marker_size );
hold on
handle_plot_vs_n(2) = loglog( n_array, Zimmermanns_times_vs_n, 's', ...
    'LineStyle', style_array(2), 'Color', color_array(2,:), 'LineWidth', linewidth, ...
    'MarkerEdgeColor', color_array(2,:), 'MarkerFaceColor', color_array(2,:), ...
    'MarkerSize', marker_size );
handle_plot_vs_n(3) = loglog( n_array, SSAF_times_vs_n, 'v', ...
    'LineStyle', style_array(3), 'Color', color_array(3,:), 'LineWidth', linewidth, ...
    'MarkerEdgeColor', color_array(3,:), 'MarkerFaceColor', color_array(3,:), ...
    'MarkerSize', marker_size );
% handle_plot_vs_p(4) = loglog( n_array(7:end), 5e-5*n_array(7:end).^(1/2), ...
%     'LineStyle', style_array(4), 'Color', gray2, 'LineWidth', linewidth );
grid on
xlabel( '$n$', 'FontSize', 16 )
ylabel( 'Avg. computational time', 'FontSize', 16 )

% Change line width of markers' edges:
drawnow;
for k=1:size(handle_plot_vs_n,2)
    handle_plot_vs_n(k).MarkerHandle.LineWidth = myMarkerLineWidth;
end

% Legend
handleLegend_plot_vs_n = legend( handle_plot_vs_n, ...
    {'[5, Alg. 1]', '[30, Alg. 1]', 'SSAF'}, ...
    'FontSize', 16, 'Location', 'NW' );

% Change line width of markers' edges in the legend:
drawnow;
for i=1:size(handle_plot_vs_n,2)
    lineEntry = findobj(handleLegend_plot_vs_n.EntryContainer, 'Object', handle_plot_vs_n(i) );
    entryMarker = findobj(lineEntry.Icon.Transform, 'Description','Icon Marker');
    entryMarker.LineWidth = myMarkerLineWidth;
end
%--------------------------------------------------------------------------
%--------------------------------------------------------------------------
%--------------------------------------------------------------------------

% Similar plot, but as a function of p:

% doubling values of p, from 2 to 256:
p_array = 2.^(1:8);

Bryners_times_vs_p = [ 0.00353, 0.00533, 0.00711, 0.01173, 0.02912, 0.08762, 0.40437, 1.94025 ];
Zimmermanns_times_vs_p = [ 0.00103, 0.00156, 0.00182, 0.00369, 0.01354, 0.03582, 0.10052, 0.47720 ];
SSAF_times_vs_p = [ 0.00086, 0.00128, 0.00115, 0.00173, 0.00453, 0.01150, 0.05657, 0.25847 ];
%--------------------------------------------------------------------------

nexttile
handle_plot_vs_p(1) = loglog( p_array, Bryners_times_vs_p, 'o', ...
    'LineStyle', style_array(1), 'Color', color_array(1,:), 'LineWidth', linewidth, ...
    'MarkerEdgeColor', color_array(1,:), 'MarkerFaceColor', color_array(1,:), ...
    'MarkerSize', marker_size );
hold on
handle_plot_vs_p(2) = loglog( p_array, Zimmermanns_times_vs_p, 's', ...
    'LineStyle', style_array(2), 'Color', color_array(2,:), 'LineWidth', linewidth, ...
    'MarkerEdgeColor', color_array(2,:), 'MarkerFaceColor', color_array(2,:), ...
    'MarkerSize', marker_size );
handle_plot_vs_p(3) = loglog( p_array, SSAF_times_vs_p, 'v', ...
    'LineStyle', style_array(3), 'Color', color_array(3,:), 'LineWidth', linewidth, ...
    'MarkerEdgeColor', color_array(3,:), 'MarkerFaceColor', color_array(3,:), ...
    'MarkerSize', marker_size );
% handle_plot_vs_p(4) = loglog( p_array(5:end), 1e-8*p_array(5:end).^3, ...
%     'LineStyle', style_array(4), 'Color', gray2, 'LineWidth', linewidth );
grid on
xlabel( '$p$', 'FontSize', 16 )
ylabel( 'Avg. computational time', 'FontSize', 16 )

% Change line width of markers' edges:
drawnow;
for k=1:size(handle_plot_vs_p,2)-1
    handle_plot_vs_p(k).MarkerHandle.LineWidth = myMarkerLineWidth;
end

% Legend
handleLegend_plot_vs_p = legend( handle_plot_vs_p, ...
    {'[5, Alg. 1]', '[30, Alg. 1]', 'SSAF'}, ...
    'FontSize', 16, 'Location', 'NW' );

% Change line width of markers' edges in the legend:
drawnow;
for i=1:size(handle_plot_vs_p,2)-1
    lineEntry = findobj(handleLegend_plot_vs_p.EntryContainer, 'Object', handle_plot_vs_p(i) );
    entryMarker = findobj(lineEntry.Icon.Transform, 'Description','Icon Marker');
    entryMarker.LineWidth = myMarkerLineWidth;
end


%--------------------------------------------------------------------------
% SAVE IMAGE TO PDF FILE
%--------------------------------------------------------------------------
fprintf('+--------------------------------------------------------------+\n');
fprintf('|                           Save plot                          |\n');
fprintf('+--------------------------------------------------------------+\n');
fileName_plot = 'plots/Comparison_avg_comput_times';
saveas( gcf, fileName_plot, 'epsc');
fprintf(' Saved graph to file %s.eps.\n', fileName_plot);
%--------------------------------------------------------------------------
