function anaPlot(data_t, data_ori, data_ind, titletext, plottype)
%ANAPLOT 绘制MIMO系统的数据图
% data_ori和data_ind的格式: {data_1, ... data_n}, 其中每个data_i @ (samples, sizes)
% 则结果一共有n行, sizes列.
    
    % 参数
    cell_size = length(data_ori);
    data_size = size(data_ori{1}, 2);
    % 检查
    if length(data_ori) ~= length(data_ind), return; end

    % 定义绘图函数
    switch plottype
        case 'plot', pfunc = @(x1, x2, x3, x4) plot(x1, x2, x3, x4);
        case 'semilogx', pfunc = @(x1, x2, x3, x4) semilogx(x1, x2, x3, x4);
        otherwise, pfunc = @(x1, x2, x3, x4) plot(x1, x2, x3, x4);
    end

    % 绘图
    figure;
    tl = tiledlayout(cell_size, data_size);
    tl.Title.String = titletext;
    % 子图
    for iter_row = 1:cell_size
        for iter_col = 1:data_size
            nexttile((iter_row-1)*data_size+iter_col);
            pfunc(data_t, data_ori{iter_row}(:, iter_col), 'Color', '#0072BD');
            orilegend = "Original System @ (" + iter_row + ", " + iter_col + ")";
            hold on;
            pfunc(data_t, data_ind{iter_row}(:, iter_col), 'Color', '#D95319');
            indlegend = "Identified System @ (" + iter_row + ", " + iter_col + ")";
            legend(orilegend, indlegend);
            hold off;
        end
    end
    
end

