function plotLegend(texts, font_size, line_width, location)
%PLOTLEGEND 图例绘制

    [lgd, lgdobj, ~, ~] = legend(texts, 'FontName', 'Arial', 'FontSize', font_size, 'FontWeight', 'bold');
    lgd.Location = location; lgd.Box = 'off';
    child_line = findobj(lgdobj, 'type', 'Line');
    set(child_line, 'LineWidth', line_width);
end

