function fig = anaPlotHistogram(cov_errors)
%ANAPLOTHISTOGRAM 绘制(对数的多重)直方图

    % 参数计算
    data_size = size(cov_errors, 1);

    % 参数定义
    legend_text = {'discrete-cISSIM (full)', 'discrete-cISSIM (reduce)', 'SIM'};
    % line_shape = {'-x', '-x', '-^'};
    line_color = {[0.8500 0.3250 0.0980], [0.4660 0.6740 0.1880], [0.9290 0.6940 0.1250]};
    font_size = 10; line_width = 1.5;

    % edges计算
    log_base = 10;
    edges_log_min = floor(log(min(cov_errors, [], 'all'))/log(log_base));
    edges_log_max = ceil(log(max(cov_errors, [], 'all'))/log(log_base));
    edges_log = edges_log_min:edges_log_max;
    edges_normal = log_base.^edges_log;
    centers_log = (edges_log(1:end-1)+edges_log(2:end))./2;
    % bins计算
    bin_counts = zeros(data_size, length(edges_log)-1);
    for iter_data = 1:data_size
        [bin_counts(iter_data, :), ~] = histcounts(cov_errors(iter_data, :), edges_normal);
    end

    % 图窗
    fig = figure;
    fig.Units = 'centimeters'; fig.Position = [0 0 14 7];
    % 绘图
    ax = axes(fig);
    hb = bar(ax, centers_log, bin_counts.');
    for iter_data = 1:data_size
        hb(iter_data).CData = repmat(line_color{iter_data}, [length(edges_log)-1 1]);
        hb(iter_data).FaceColor = 'flat';
        hb(iter_data).EdgeColor = [.3,.3,.3];
        hb(iter_data).BarWidth = 1;
        hb(iter_data).LineWidth = line_width;
    end
    % x轴
    edges_txt = repmat("10^{", size(edges_log)) + string(edges_log) + repmat("}", size(edges_log));
    xticks(ax, edges_log);
    xticklabels(ax, edges_txt);

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

    % 坐标轴修改
    xlim(ax, [edges_log(1) edges_log(end)]);
    % 显示标签
    anaPlotLabel('\eta_{s}', 'Frequency', font_size);
    % 显示图例
    anaPlotLegend(legend_text(1:data_size), font_size, line_width, 'northwest');

end

% [~, edges_log] = histcounts(log10(sum(cov_errors, 1)));
% edges_normal = 10.^edges_log;
% histogram(ax, cov_errors(iter_data, :), edges_normal, 'EdgeColor', line_color{iter_data}, 'FaceColor', line_color{iter_data});
% set(ax, 'xscale','log');

