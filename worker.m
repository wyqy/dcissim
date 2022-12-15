%% SISO
%#ok<*NASGU> 
clear
save_valid = true; para_experiment_count = 100;
main_dcissim_siso

%% MIMO - recursive
%#ok<*NASGU> 
clear
save_valid = true;
main_dcissim_recursive


%% MIMO - covariance
% 需要注释掉idenCovariance.m中的if min(eig(cov_zr)) < 0 || min(eig(cov_zrt_2)) < 0判断句
%#ok<*NASGU> 
clear
save_valid = true; para_experiment_count = 100;
main_dcissim_als
