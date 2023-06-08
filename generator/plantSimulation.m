function ret_struct = plantSimulation(varargin)
%PLANTSIMULATION 返回(开环)系统带噪声的仿真结果

    % 输入解析
    parser = inputParser;
    % 输入定义
    addParameter(parser, 'type', 'fixed', @(i)(ischar(i)));
    addParameter(parser, 'x_size', 1, @(i)(isnumeric(i)&&isscalar(i)));
    addParameter(parser, 'y_size', 1, @(i)(isnumeric(i)&&isscalar(i)));
    addParameter(parser, 'u_size', 1, @(i)(isnumeric(i)&&isscalar(i)));
    addParameter(parser, 'samples_period', 1, @(i)(isnumeric(i)&&isscalar(i)));
    addParameter(parser, 'excitation', 1, @(i)(isnumeric(i)&&ismatrix(i)));
    addParameter(parser, 'snr', nan, @(i)(isnumeric(i)&&isvector(i)));
    % 输入提取
    parse(parser, varargin{:});
    type = parser.Results.type;  % 系统参数来源
    x_size = parser.Results.x_size;  % x size
    y_size = parser.Results.y_size;  % y size
    u_size = parser.Results.u_size;  % u size
    samples_period = parser.Results.samples_period;  % 单周期采样点数
    excitation = parser.Results.excitation;  % 激励信号
    xyu_snr = parser.Results.snr;  % x,y,u聚合的信噪比

    % 参数构造
    switch type
        case 'fixed'
            A = [0.8 0; 0 0.2];
            B = [1; 1];
            C = [1 0; 1 1];
            D = [0; 0];
            seed = rng().Seed;
            covariance = blkdiag([1 0; 0 1], [1 0; 0 1], 1);
            plant_info = struct('A', A, 'B', B, 'C', C, 'D', D, ...
                'seed', seed, 'cov', covariance);
        case 'generated'
            % 生成随机种子
            seed = rng().Seed;
            rs = RandStream('dsfmt19937', 'Seed', seed);
            % 随机生成系统参数
            is_feasible = false;
            while ~is_feasible
                rand_sys = drss(x_size, y_size, u_size);
                is_feasible = max(abs(eig(rand_sys.A, 'vector'))) <= 0.95;  % 经验值
                if y_size >= x_size  % 用于方差估计
                    is_cfullrank = rank(rand_sys.C) >= x_size;
                    is_feasible = is_feasible && is_cfullrank;
                end
            end
            % 随机生成协方差矩阵
            scale_cov = 5;
            covariance_x = scale_cov*semidefMatrixBuilder(rs, x_size);
            covariance_y = scale_cov*semidefMatrixBuilder(rs, y_size);
            covariance_u = scale_cov*semidefMatrixBuilder(rs, u_size);
            covariance = blkdiag(covariance_x, covariance_y, covariance_u);
            % 构造struct
            plant_info = struct('A', rand_sys.A, 'B', rand_sys.B, 'C', rand_sys.C, 'D', rand_sys.D, ...
                'seed', seed, 'cov', covariance);
        otherwise, plant_info = load('parameter\paraPlant.mat', 'para'); plant_info = plant_info.para; seed = plant_info.seed;
    end
    % 参数计算
    x_size = size(plant_info.A, 1);
    xyu_size = x_size+size(plant_info.C, 1)+size(plant_info.B, 2);
    signal_sample = size(excitation, 2);

    % 试生成仿真数据
    x_init_powertest = zeros(x_size, 1);
    uk_powertest = excitation(:, 1:samples_period);
    noise_powertest = zeros(xyu_size, samples_period);
    [xk_test, yk_test, uk_test] = plantModel(plant_info, x_init_powertest, uk_powertest, noise_powertest);
    xyuk_test = [xk_test; yk_test; uk_test];
    xyuk_test = xyuk_test(:, end-samples_period+1:end);
    % 修正协方差矩阵
    covariance = genCovariance(plant_info, xyu_snr, xyuk_test);
    nozero_locs = diag(covariance) ~= 0;
    scale_transfer = zeros(xyu_size, xyu_size);
    scale_transfer(nozero_locs, nozero_locs) = ctranspose(chol(covariance(nozero_locs, nozero_locs)));
    
    % 计算噪声信号
    rs = RandStream.create('mrg32k3a', 'NumStreams', xyu_size, 'Seed', seed, 'CellOutput', true);
    noise = zeros(xyu_size, signal_sample);
    for iter = 1:xyu_size, noise(iter, :) = randn(rs{iter}, 1, signal_sample); end
    noise = scale_transfer*noise;
    % 计算等效Kalman增益
    moment_xx = dlyap(plant_info.A, covariance(1:x_size, 1:x_size));
    moment_yy = plant_info.C*moment_xx*plant_info.C.' + covariance(x_size+1:x_size+y_size, x_size+1:x_size+y_size);
    moment_xy = plant_info.A*moment_xx*plant_info.C.' + covariance(1:x_size, x_size+1:x_size+y_size);
    [~, kalman_gain, ~] = idare(plant_info.A.', -plant_info.C.', [], -moment_yy, moment_xy, []);
    kalman_gain = -kalman_gain.';

    % 更新方差数据
    plant_info.cov = covariance;
    plant_info.covariance_struct = struct('cov_all', covariance, 'kalman', kalman_gain, ...
        'cov_xx', covariance(1:x_size, 1:x_size), ...
        'cov_xy', covariance(1:x_size, x_size+1:x_size+y_size), ...
        'cov_yy', covariance(x_size+1:x_size+y_size, x_size+1:x_size+y_size), ...
        'cov_xu', covariance(1:x_size, x_size+y_size+1:end), ...
        'cov_yu', covariance(x_size+1:x_size+y_size, x_size+y_size+1:end), ...
        'cov_uu', covariance(x_size+y_size+1:end, x_size+y_size+1:end));

    % 生成仿真数据
    x_init = zeros(x_size, 1);  % 零初始化
    [xk, yk, uk] = plantModel(plant_info, x_init, excitation, noise);
    % 返回值
    ret_struct = struct('plant_info', plant_info, 'xk', xk, 'yk', yk, 'uk', uk, 'noise', noise, 'period_samples', samples_period);

end


function ret_mat = semidefMatrixBuilder(rand_rs, n)
% 正定矩阵生成(0~1)
    ret_mat = rand(rand_rs, n, n);
    ret_mat = 0.5.*(ret_mat+ret_mat.');
    ret_mat = ret_mat + n.*eye(n);
end

