function ret_struct = genPlantSimulation(varargin)
%GENPLANTSIMULATION 返回(开环)系统带噪声的仿真结果

    % 输入解析
    parser = inputParser;
    % 输入定义
    addParameter(parser, 'period_samples', 1, @(i)(isnumeric(i)&&isscalar(i)));
    addParameter(parser, 'excitation', 1, @(i)(isnumeric(i)&&ismatrix(i)));
    addParameter(parser, 'snr', nan, @(i)(isnumeric(i)&&isvector(i)));
    % 输入提取
    parse(parser, varargin{:});
    period_samples = parser.Results.period_samples;  % 单周期采样点数
    excitation = parser.Results.excitation;  % 激励信号
    xyu_snr = parser.Results.snr;  % x,y,u聚合的信噪比

    % 参数提取
    plant_info = load('parameter\paraPlant.mat', 'para');
    plant_info = plant_info.para;
    
    % 参数计算
    x_size = size(plant_info.A, 1);
    xyu_size = x_size+size(plant_info.C, 1)+size(plant_info.B, 2);
    samples = size(excitation, 2);

    % 试生成仿真数据
    repeat_count = 3;  % 自定义参数
    x_init_powertest = zeros(x_size, 1);
    un_powertest = excitation(:, 1:repeat_count*period_samples);
    noise_powertest = zeros(xyu_size, repeat_count*period_samples);
    [xn_test, yn_test, un_test] = genPlantModel(plant_info, x_init_powertest, un_powertest, noise_powertest);
    % 生成噪声信号
    xyun_test = [xn_test; yn_test; un_test];
    xyun_test = xyun_test(:, end-period_samples+1:end);
    noise = genNoiser(plant_info, samples, xyu_snr, xyun_test);

    % 生成仿真数据
    x_init = zeros(x_size, 1);
    [xn, yn, un] = genPlantModel(plant_info, x_init, excitation, noise);
    
    % 返回值
    ret_struct = struct('xn', xn, 'yn', yn, 'un', un, 'noise', noise, 'period_samples', period_samples);

end

