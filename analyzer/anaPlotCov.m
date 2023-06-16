function fig = anaPlotCov(cov_norm)
%ANASEMILOGX 绘制方差范数对比图 - 使用半对数坐标系

    % 参数计算
    data_size = size(cov_norm, 1);
    cov_xdim = 0:size(cov_norm, 2)-1;

    % 参数定义
    legend_text = {'Original', 'discrete-cISSIM (full)', 'discrete-cISSIM (reduce)', 'SIM'};
    line_shape = {'-.o', '-x', '-x', '-^'};
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
    % axis(ax, [1 size(cov_norm, 2) -inf inf]);

    % 寻找axes
    child_axes = findobj(fig.Children, 'type', 'Axes');
    for iter_axes = 1:length(child_axes)
        % 优化格式
        plotAux(child_axes(iter_axes), font_size, line_width);
        % 寻找line
        child_line = findobj(child_axes(iter_axes), 'type', 'Line');
        for iter_line = 1:length(child_line)
            plotAux(child_line(iter_line), font_size, line_width);
        end
    end

    % 显示标签
    plotLabel('Lag of \color[rgb]{1,1,1} o \color[rgb]{.3,.3,.3} {\Xi}_{rr}^{i}', 'H_{2} Norm', font_size);
    % 显示图例
    plotLegend(legend_text(1:data_size), font_size, line_width, 'north');

    % 保存图片
    % filepath = 'fig';
    % saveas(gcf, [filepath 'steering_lateral_bode_plot.eps'], 'epsc')

end

