function fig = plotBode(ss_cell, varargin)
%PLOTBODE 绘制Bode图
    
    % 输入定义
    parser = inputParser;
    addParameter(parser, 'sample', 1, @(i)(isnumeric(i)&&isscalar(i)));
    parse(parser, varargin{:});
    sample_step = parser.Results.sample;  % 采样速率
    % 参数计算
    model_size = length(ss_cell);
    y_size = size(ss_cell{1}.C, 1);
    u_size = size(ss_cell{1}.B, 2);
    % 兼容自定义格式
    for iter_model = 1:model_size
        if ~(isa(ss_cell{iter_model}, 'ss') || isa(ss_cell{iter_model}, 'idss'))
            ss_cell{iter_model} = ss(ss_cell{iter_model}.A, ss_cell{iter_model}.B, ss_cell{iter_model}.C, ss_cell{iter_model}.D, sample_step);
        end
    end

    % 参数定义
    legend_text = {'Original', 'discrete-cISSIM (full)', 'discrete-cISSIM (reduce)', 'SIM'};
    line_shape = {'-.', '-', '-', '-'};
    line_color = {'#0072BD', '#D95319', '#77AC30', '#EDB120'};
    font_size = 10; line_width = 1.5;

    % 图窗
    fig = figure;
    fig.Units = 'centimeters'; fig.Position = [0 0 14 10];
    
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
        plotoptions.XLabel.String = '\fontname{Arial}{Frequency}'; plotoptions.XLabel.FontSize = font_size; plotoptions.XLabel.FontWeight = 'bold'; plotoptions.XLabel.Interpreter = 'tex';
        plotoptions.YLabel.String = {'\fontname{Arial}{Amplitude}', '\fontname{Arial}{Phase}'}; plotoptions.YLabel.FontSize = font_size; plotoptions.YLabel.FontWeight = 'bold'; plotoptions.YLabel.Interpreter = 'tex';
        if y_size ~= 1 || u_size ~= 1, plotoptions.Title.String = ['\fontname{Arial}u' int2str(location_row) ' - y' int2str(locaiton_col)]; plotoptions.Title.FontSize = font_size;
        else, plotoptions.Title.String = ''; end
        plotoptions.Title.FontSize = font_size; plotoptions.Title.FontWeight = 'bold'; plotoptions.Title.Interpreter = 'tex';
        temp_state = state_cell_row; temp_state{location_row} = 'on'; plotoptions.InputVisible = temp_state;
        temp_state = state_cell_col; temp_state{locaiton_col} = 'on'; plotoptions.OutputVisible = temp_state;
        % 绘制Bode图
        % bodeplot(ax, ss_original, 'b-', ss_identified, 'r-', {0.1, 10}, plotoptions);
        % bodeplot(ax, ss_original, 'b-', ss_identified, 'r-', plotoptions);
        for iter_model = 1:model_size
            bodeplot(ax, ss_cell{iter_model}, line_shape{iter_model}, plotoptions); hold on;
        end
    end

    % 删除标题
    if y_size == 1 && u_size == 1
        child_axes = findobj(fig.Children, 'type', 'Axes');
        for iter_axes = 1:length(child_axes)
            delete(child_axes(iter_axes).Title);
        end
    end
    
    % 寻找axes
    child_axes = findobj(fig.Children, 'type', 'Axes');
    for iter_axes = 1:length(child_axes)
        % 优化格式
        plotAux(child_axes(iter_axes), font_size, line_width);
        % 寻找line
        child_line = findobj(child_axes(iter_axes), 'type', 'Line');
        for iter_line = 1:length(child_line)
            % 优化格式
            plotAux(child_line(iter_line), font_size, line_width);
        end
        % 修改颜色(i dont know why)
        line_location = 1;
        if length(child_line)-1 == y_size*model_size
            for iter_line = length(child_line):-y_size:2
                for iter_dy = 0:y_size-1, child_line(iter_line-iter_dy).Color = line_color{line_location}; end
                line_location = line_location + 1;
            end
        elseif length(child_line)-1 == u_size*model_size
            for iter_line = length(child_line):-u_size:2
                for iter_dy = 0:u_size-1, child_line(iter_line-iter_dy).Color = line_color{line_location}; end
                line_location = line_location + 1;
            end
        end
        
    end

    % 在最后一个图上显示图例
    plotLegend(legend_text(1:model_size), font_size, line_width, 'southwest');
    
    % 保存图片
    % filepath = 'fig';
    % saveas(gcf, [filepath 'steering_lateral_bode_plot.eps'], 'epsc')
end



