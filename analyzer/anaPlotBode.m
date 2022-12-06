function fig = anaPlotBode(ss_original, ss_identified, varargin)
%ANABODE 绘制Bode图
    
    % 输入定义
    parser = inputParser;
    addParameter(parser, 'sample', 1, @(i)(isnumeric(i)&&isscalar(i)));
    parse(parser, varargin{:});
    sample_step = parser.Results.sample;  % 采样速率
    % 兼容自定义格式
    if ~(isa(ss_original, 'ss') || isa(ss_original, 'idss'))
        ss_original = ss(ss_original.A, ss_original.B, ss_original.C, ss_original.D, sample_step);
    end
    if ~(isa(ss_identified, 'ss') || isa(ss_identified, 'idss'))
        ss_identified = ss(ss_identified.A, ss_identified.B, ss_identified.C, ss_identified.D, sample_step);
    end

    % 参数定义
    orilegend = 'Original';
    indlegend = 'Identified';
    font_size = 10;
    y_size = size(ss_original.C, 1);
    u_size = size(ss_original.B, 2);

    % 图窗
    fig = figure;
    % f = figure; f.Units = 'centimeters'; f.Position = [0 0 14 10];
    
    % 变量预定义
    state_cell_row = cell(u_size, 1);  % 行表示输入
    state_cell_col = cell(y_size, 1);  % 列表示输出
    for iter_init = 1:u_size, state_cell_row{iter_init} = 'off'; end
    for iter_init = 1:y_size, state_cell_col{iter_init} = 'off'; end
    
    for iter_ax = 1:y_size*u_size
        % 子图位置
        ax = subplot(u_size, y_size, iter_ax);
        % 参数计算
        location_row = ceil(iter_ax/y_size);  % 同时是输入标号
        locaiton_col = mod(iter_ax-1, y_size)+1;  % 同时是输出标号

        % Bode图配置
        plotoptions = bodeoptions('cstprefs');  % ctrlpref
        plotoptions.PhaseWrapping = 'on';
        plotoptions.IOGrouping = 'all';
        plotoptions.XLabel.String = '\fontname{Cambria}{Frequency}'; plotoptions.XLabel.FontSize = font_size; plotoptions.XLabel.FontWeight = 'bold'; plotoptions.XLabel.Interpreter = 'tex';
        plotoptions.YLabel.String = {'\fontname{Cambria}{Amplitude}', '\fontname{Cambria}{Phase}'}; plotoptions.YLabel.FontSize = font_size; plotoptions.YLabel.FontWeight = 'bold'; plotoptions.YLabel.Interpreter = 'tex';
        plotoptions.Title.String = ['\fontname{Cambria}u' int2str(location_row) ' - y' int2str(locaiton_col)]; plotoptions.Title.FontSize = font_size; plotoptions.Title.FontWeight = 'bold'; plotoptions.Title.Interpreter = 'tex';
        temp_state = state_cell_row; temp_state{location_row} = 'on'; plotoptions.InputVisible = temp_state;
        temp_state = state_cell_col; temp_state{locaiton_col} = 'on'; plotoptions.OutputVisible = temp_state;
        % 绘制Bode图
        % bodeplot(ax, ss_original, 'b-', ss_identified, 'r-', {0.1, 10}, plotoptions);
        bodeplot(ax, ss_original, 'b-', ss_identified, 'r-', plotoptions);
        % 优化格式
        anaPlotAux(fig.Children(2), font_size); anaPlotAux(fig.Children(3), font_size);
        % 在第一个图上显示图例
        if iter_ax == 1, bodeLabelCaller(orilegend, indlegend, font_size); end
    end
    
    % 保存图片
    % filepath = 'fig';
    % saveas(gcf, [filepath 'steering_lateral_bode_plot.eps'], 'epsc')
end

function bodeLabelCaller(text1, text2, font_size)
    lgd = legend(text1, text2, 'FontName', 'Cambria', 'FontSize', font_size, 'FontWeight', 'bold');
    lgd.Location = 'northwest'; lgd.Box = 'off';
end


