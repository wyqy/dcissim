function anaPlotLabel(xtexts, ytexts, font_size)
%ANAPLOTLABEL 标签绘制

    xlabel(xtexts, 'FontName', 'Cambria', 'FontSize', font_size, 'FontWeight', 'bold');
    ylabel(ytexts, 'FontName', 'Cambria', 'FontSize', font_size, 'FontWeight', 'bold');
end
