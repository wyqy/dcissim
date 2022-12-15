function fig = anaPlotRecursive(cov_norm, norm_period)
%ANAPLOTRECURSIVE 绘制迭代的范数变换图, 没有标记点

    % 参数计算
    data_size = size(cov_norm, 1);
    cov_xdim = (1:size(cov_norm, 2))./norm_period;

    % 参数定义
    legend_text = {'Original', 'discrete-cISSIM (full)', 'discrete-cISSIM (reduce)', 'SIM'};
    line_shape = {'-.', '-', '-', '-'};
    line_color = {'#0072BD', '#D95319', '#77AC30', '#EDB120'};
    font_size = 10; line_width = 1.5;

    % 图窗
    fig = figure;
    fig.Units = 'centimeters'; fig.Position = [0 0 14 7];
    % 绘图
    ax = axes(fig);
    for iter_data = 1:data_size
        plot(ax, cov_xdim, cov_norm(iter_data, :), line_shape{iter_data}, 'Color', line_color{iter_data}); hold on;
    end
    xlabel(ax, '')

    % 寻找axes
    child_axes = findobj(fig.Children, 'type', 'Axes');
    for iter_axes = 1:length(child_axes)
        % 优化格式
        anaPlotAux(child_axes(iter_axes), font_size, line_width);
        % 寻找line
        child_line = findobj(child_axes(iter_axes), 'type', 'Line');
        for iter_line = 1:length(child_line)
            anaPlotAux(child_line(iter_line), font_size, line_width);
        end
    end

    % 显示标签
    anaPlotLabel('Period', 'H_{2} Norm', font_size);
    % 显示图例
    anaPlotLegend(legend_text(1:data_size), font_size, line_width, 'north');

end

