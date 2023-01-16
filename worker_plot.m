%% SISO load
clear variables
load('results\mat_auto_siso.mat');

%% SISO Bode plot
analysis_bode_location = 48;
analysis_bode_cell = {result_original_cell{analysis_bode_location}, result_full_cell{analysis_bode_location}, result_reduce_cell{analysis_bode_location}, result_sim_cell{analysis_bode_location}};
fig_bode = anaPlotBode(analysis_bode_cell, 'sample', para_sim_step); % sgtitle(fig, 'Bode plot for minimum d\_{Hinf}');
exportgraphics(fig_bode, '..\..\%latex%\derivation - mod10\siso_bode_1.eps')

%% SISO covariance plot
analysis_cov_location = 48;
fig_cov = anaPlotCov(analysis_outcov_series(:, :, analysis_cov_location));
exportgraphics(fig_cov, '..\..\%latex%\derivation - mod10\siso_covariance_1.eps')

%% SISO covariance histogram
% fig_cov = anaPlotCov(analysis_outcov_series(:, :, analysis_cov_location));
fig_hist = anaPlotHistogram(analysis_outcov_error_norm);
exportgraphics(fig_hist, '..\..\%latex%\derivation - mod10\siso_covariance_hist_1.eps')

%% MIMO recursive load
clear variables
load('results\mat_auto_recursive.mat');

%% MIMO recursive plot
fig_recursive = anaPlotRecursive(analysis_h2, 10);
exportgraphics(fig_recursive, '..\..\%latex%\derivation - mod10\mimo_recursive_1.eps')

%% MIMO scatter load
clear variables
load('results\mat_auto_als.mat');

%% MIMO scatter plot
fig_scatter = anaPlotScatter(analysis_spot_inout_3d);
set(fig_scatter,'renderer','Painters');
exportgraphics(fig_scatter, '..\..\%latex%\derivation - mod10\mimo_covariance_1.eps')


