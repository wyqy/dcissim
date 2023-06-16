function plotLabel(xtexts, ytexts, font_size)
%PLOTLABEL 标签绘制

    xlabel(xtexts, 'FontName', 'Arial', 'FontSize', font_size, 'FontWeight', 'bold');
    ylabel(ytexts, 'FontName', 'Arial', 'FontSize', font_size, 'FontWeight', 'bold');
end
