function ret_struct = plantSimulation(varargin)
%PLANTSIMULATION 返回(开环)系统带噪声的仿真结果

    % 输入解析
    parser = inputParser;
    % 输入定义
    addParameter(parser, 'type', 'fixed', @(i)(ischar(i)));
    addParameter(parser, 'period_samples', 1, @(i)(isnumeric(i)&&isscalar(i)));
    addParameter(parser, 'excitation', 1, @(i)(isnumeric(i)&&ismatrix(i)));
    addParameter(parser, 'snr', nan, @(i)(isnumeric(i)&&isvector(i)));
    addParameter(parser, 'x_size', 1, @(i)(isnumeric(i)&&isscalar(i)));
    addParameter(parser, 'y_size', 1, @(i)(isnumeric(i)&&isscalar(i)));
    addParameter(parser, 'u_size', 1, @(i)(isnumeric(i)&&isscalar(i)));
    % 输入提取
    parse(parser, varargin{:});
    type = parser.Results.type;  % 系统参数来源
    period_sample = parser.Results.period_samples;  % 单周期采样点数
    excitation = parser.Results.excitation;  % 激励信号
    xyu_snr = parser.Results.snr;  % x,y,u聚合的信噪比
    x_size = parser.Results.x_size;  % x size
    y_size = parser.Results.y_size;  % y size
    u_size = parser.Results.u_size;  % u size

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
        otherwise, plant_info = load('parameter\paraPlant.mat', 'para'); plant_info = plant_info.para;
    end
    % 参数计算
    x_size = size(plant_info.A, 1);
    xyu_size = x_size+size(plant_info.C, 1)+size(plant_info.B, 2);
    signal_sample = size(excitation, 2);

    % 试生成仿真数据
    x_init_powertest = zeros(x_size, 1);
    un_powertest = excitation(:, period_sample+1:2*period_sample);
    noise_powertest = zeros(xyu_size, period_sample);
    [xn_test, yn_test, un_test] = plantModel(plant_info, x_init_powertest, un_powertest, noise_powertest);
    % 生成噪声信号
    xyun_test = [xn_test; yn_test; un_test];
    xyun_test = xyun_test(:, end-period_sample+1:end);
    [noise, covariance] = genNoiser(plant_info, signal_sample, xyu_snr, xyun_test);
    % 更新方差数据
    plant_info.cov = covariance;

    % 生成仿真数据
    x_init = zeros(x_size, 1);  % 零初始化
    [xn, yn, un] = plantModel(plant_info, x_init, excitation, noise);
    % 返回值
    ret_struct = struct('plant_info', plant_info, 'xn', xn, 'yn', yn, 'un', un, 'noise', noise, 'period_samples', period_sample);

end


function ret_mat = semidefMatrixBuilder(rand_rs, n)
% 正定矩阵生成(0~1)
    ret_mat = rand(rand_rs, n, n);
    ret_mat = 0.5.*(ret_mat+ret_mat.');
    ret_mat = ret_mat + n.*eye(n);
end

