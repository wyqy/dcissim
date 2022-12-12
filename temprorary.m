%% SISO prepare
load('results\mat_full_auto_siso.mat', ...
    'analysis_error_h2', 'analysis_error_hinf', 'result_original_cell', 'result_dcissim_cell', 'result_sim_cell', ...
    'analysis_outcov_error_norm', 'analysis_outcov_series', 'para_sim_step');
analysis_error_h2_all = analysis_error_h2;
analysis_error_hinf_all = analysis_error_hinf;
analysis_outcov_error_norm_all = analysis_outcov_error_norm;
analysis_outcov_series_all = analysis_outcov_series;
result_full_dcissim_cell = result_dcissim_cell;

load('results\mat_reduce_auto_siso.mat', ...
    'analysis_error_h2', 'analysis_error_hinf', 'result_dcissim_cell', ...
    'analysis_outcov_series');
analysis_error_h2_all = [analysis_error_h2_all(1, :); analysis_error_h2(1, :); analysis_error_h2_all(2, :)];
analysis_error_hinf_all = [analysis_error_hinf_all(1, :); analysis_error_hinf(1, :); analysis_error_hinf_all(2, :)];
analysis_outcov_error_norm_all = [analysis_outcov_error_norm_all(1, :); analysis_outcov_error_norm(1, :); analysis_outcov_error_norm_all(2, :)];
analysis_outcov_series_all = [analysis_outcov_series_all(1:2, :, :); analysis_outcov_series(2, :, :); analysis_outcov_series_all(3, :, :)];
result_reduce_dcissim_cell = result_dcissim_cell;

result_reduce_dcissim_cell = sparsedMat(result_reduce_dcissim_cell);
result_full_dcissim_cell = sparsedMat(result_full_dcissim_cell);

save('results\mat_auto_siso_processed.mat', '-v7.3', ...
    'analysis_error_h2_all', 'analysis_error_hinf_all', 'result_original_cell', 'result_full_dcissim_cell', 'result_reduce_dcissim_cell', 'result_sim_cell', ...
    'analysis_outcov_error_norm_all', 'analysis_outcov_series_all', 'para_sim_step');

%% SISO load 
load('results\mat_auto_siso_processed.mat');

%% SISO Bode plot
analysis_bode_location = 48;
analysis_bode_cell = {result_original_cell{analysis_bode_location}, result_full_dcissim_cell{analysis_bode_location}, result_reduce_dcissim_cell{analysis_bode_location}, result_sim_cell{analysis_bode_location}};
fig = anaPlotBode(analysis_bode_cell, 'sample', para_sim_step); % sgtitle(fig, 'Bode plot for minimum d\_{Hinf}');
exportgraphics(gcf, '..\..\%latex%\derivation - mod3\fig\siso_bode_1.eps')

%% SISO covariance plot
analysis_cov_location = 48;
ax = anaPlotCov(analysis_outcov_series_all(:, :, analysis_cov_location));
exportgraphics(gcf, '..\..\%latex%\derivation - mod3\fig\siso_covariance_1.eps')


%% aux function
function in_cell = sparsedMat(in_cell)
    cell_size = length(in_cell);
    for iter_cell = 1:cell_size
        in_cell{iter_cell}.S = sparse(in_cell{iter_cell}.S);
    end
end