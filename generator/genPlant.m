function plant_info = genPlant(varargin)
%GENPLANT 构造系统

    % 输入解析
    parser = inputParser;
    % 输入定义
    addParameter(parser, 'type', 'fixed', @(i)(ischar(i)));
    addParameter(parser, 'excitation', 1, @(i)(isnumeric(i)&&ismatrix(i)));
    addParameter(parser, 'x_size', 1, @(i)(isnumeric(i)&&isscalar(i)));
    addParameter(parser, 'y_size', 1, @(i)(isnumeric(i)&&isscalar(i)));
    addParameter(parser, 'u_size', 1, @(i)(isnumeric(i)&&isscalar(i)));
    addParameter(parser, 'samples_period', 1, @(i)(isnumeric(i)&&isscalar(i)));
    addParameter(parser, 'snr', nan, @(i)(isnumeric(i)&&isvector(i)));
    % 输入提取
    parse(parser, varargin{:});
    type = parser.Results.type;  % 系统参数来源
    excitation = parser.Results.excitation;  % 激励信号
    x_size = parser.Results.x_size;  % x size
    y_size = parser.Results.y_size;  % y size
    u_size = parser.Results.u_size;  % u size
    T = parser.Results.samples_period;  % 单周期采样点数
    snr = parser.Results.snr;  % x,y,u聚合的信噪比

    % 参数计算
    xyu_size = x_size + y_size + u_size;

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
        case 'random'
            % 生成随机种子
            seed = rng().Seed;
            rs = RandStream('dsfmt19937', 'Seed', seed);
            % 随机生成系统参数
            is_feasible = false;
            while ~is_feasible
                rand_sys = drss(x_size, y_size, u_size);
                is_feasible = max(abs(eig(rand_sys.A, 'vector'))) <= 0.8 && ...
                              min(abs(eig(rand_sys.A, 'vector'))) >= 0.2;  % 不可过大也不可过小
                is_feasible = is_feasible && rank(obsv(rand_sys)) == x_size && ...
                                             rank(ctrb(rand_sys)) == x_size; % 可控可观
                is_feasible = is_feasible && max(abs(zpk(rand_sys).Z{1})) < 1;  % 最小相位系统
                if y_size >= x_size  % 用于方差估计
                    is_cfullrank = rank(rand_sys.C) >= x_size;
                    is_feasible = is_feasible && is_cfullrank;
                end
                rand_sys.D = zeros(y_size, u_size);  % 非直通系统
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
        otherwise, plant_info = load('parameter\paraPlant.mat', 'para'); plant_info = plant_info.para;
    end

    % 试生成仿真数据
    [xk_test, yk_test, uk_test] = plantModel(plant_info, zeros(x_size, 1), excitation(:, T+1:2*T), zeros(xyu_size, T));
    xyuk_test = [xk_test; yk_test; uk_test];
    xyuk_test = xyuk_test(:, end-T+1:end);
    % 根据信噪比计算噪声功率
    signal_power = transpose(bandpower(xyuk_test.'));
    noise_power = signal_power./(10.^(snr./10));
    
    % 根据噪声功率更新方差数据
    scale_matrix = diag(sqrt(noise_power./diag(covariance)));
    covariance = scale_matrix*covariance*scale_matrix;

    % 计算不含NaN的snr的正态分布
    valid_locs = ~isnan(diag(covariance));
    temp_cov = covariance;
    covariance = zeros(xyu_size);
    covariance(valid_locs, valid_locs) = temp_cov(valid_locs, valid_locs);

    % 更新方差数据 & 构造返回值
    plant_info.cov = covariance;
    plant_info.covariance_struct = struct('cov_all', covariance, ...
        'cov_xx', covariance(1:x_size, 1:x_size), ...
        'cov_xy', covariance(1:x_size, x_size+1:x_size+y_size), ...
        'cov_yy', covariance(x_size+1:x_size+y_size, x_size+1:x_size+y_size), ...
        'cov_xu', covariance(1:x_size, x_size+y_size+1:end), ...
        'cov_yu', covariance(x_size+1:x_size+y_size, x_size+y_size+1:end), ...
        'cov_uu', covariance(x_size+y_size+1:end, x_size+y_size+1:end));

end

function ret_mat = semidefMatrixBuilder(rand_rs, n)
% 正定矩阵生成(0~1)
    ret_mat = rand(rand_rs, n, n);
    ret_mat = 0.5.*(ret_mat+ret_mat.');
    ret_mat = ret_mat + n.*eye(n);
end


