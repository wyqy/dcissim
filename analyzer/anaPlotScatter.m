function fig = anaPlotScatter(source_data)
%ANAPLOTSCATTER 绘制2d / 3d 误差分布图
% 数据维度(data_dim, 2/3, sample_number)

    % 参数计算
    data_size = size(source_data, 1);
    dim_size = size(source_data, 2);

    % 参数定义
    legend_text = {'ALS (simplified)', 'ALS (classical)'};
    spot_mark = {'x', 'x', 'x'};
    spot_size = [36 36 36];
    spot_color = {[0.8500 0.3250 0.0980], [0.4660 0.6740 0.1880], [0.9290 0.6940 0.1250]};
    font_size = 10; line_width = 1.5;

    % 图窗
    fig = figure;
    fig.Units = 'centimeters'; fig.Position = [0 0 14 10];
    % 绘图
    ax = axes(fig);
    if dim_size == 2
        for iter_data = 1:data_size
            scatter(ax, shiftdim(source_data(iter_data, 1, :)), shiftdim(source_data(iter_data, 2, :)), ...
                spot_size(iter_data), spot_color{iter_data}, spot_mark{iter_data}, 'LineWidth', line_width); hold on;
        end
        xlabel(ax, 'error of w', 'FontSize', font_size, 'FontWeight', 'bold', 'FontName', 'Arial');
        ylabel(ax, 'error of v', 'FontSize', font_size, 'FontWeight', 'bold', 'FontName', 'Arial');
    elseif dim_size == 3
        for iter_data = 1:data_size
            scatter3(ax, shiftdim(source_data(iter_data, 1, :)), shiftdim(source_data(iter_data, 2, :)), shiftdim(source_data(iter_data, 3, :)), ...
                spot_size(iter_data), spot_color{iter_data}, spot_mark{iter_data}, 'LineWidth', line_width); hold on;
        end
        xlabel(ax, 'error of \color[rgb]{1,1,1} \ldots \color[rgb]{.3,.3,.3} \omega', 'FontSize', font_size, 'FontWeight', 'bold', 'FontName', 'Arial');
        ylabel(ax, 'error of \color[rgb]{1,1,1} \ldots \color[rgb]{.3,.3,.3} \nu', 'FontSize', font_size, 'FontWeight', 'bold', 'FontName', 'Arial');
        zlabel(ax, 'error of \color[rgb]{1,1,1} \ldots \color[rgb]{.3,.3,.3} \tau', 'FontSize', font_size, 'FontWeight', 'bold', 'FontName', 'Arial');
    end

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

    % 显示图例
    anaPlotLegend(legend_text(1:data_size), font_size, line_width, 'north');

    % 保存图片
    % filepath = 'fig';
    % saveas(gcf, [filepath 'steering_lateral_bode_plot.eps'], 'epsc')

end

