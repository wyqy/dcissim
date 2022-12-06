function ax = anaPlotCov(cov_norm_original, cov_norm_indentified)
%ANASEMILOGX 绘制方差范数对比图 - 使用半对数坐标系

    % 参数定义
    orilegend = 'Original';
    indlegend = 'Identified';
    font_size = 10;

    % 图窗
    fig = figure;
    % f = figure; f.Units = 'centimeters'; f.Position = [0 0 14 10];
    % 绘图
    ax = axes(fig);
    plot(ax, cov_norm_original, '-.o'); hold on;
    plot(ax, cov_norm_indentified, '-.x'); hold on;
    anaPlotAux(ax, font_size);
    lgd = legend(orilegend, indlegend, 'FontName', 'Cambria', 'FontSize', font_size, 'FontWeight', 'bold');
    lgd.Location = 'northwest'; lgd.Box = 'off';

    % 保存图片
    % filepath = 'fig';
    % saveas(gcf, [filepath 'steering_lateral_bode_plot.eps'], 'epsc')

end

