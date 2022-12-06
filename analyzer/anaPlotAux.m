function anaPlotAux(ax, font_size)
%ANAPLOTAUX 修改axis的格式

    % 常规美化
    ax.LineWidth = 1;
    ax.Box = 'on';
    ax.TickDir = 'in';
    ax.XMinorTick = 'on';
    ax.YMinorTick = 'on';
    ax.XGrid = 'on';
    ax.YGrid = 'on';
    ax.GridLineStyle = '--';
    ax.XColor = [.3,.3,.3];
    ax.YColor = [.3,.3,.3];
    ax.FontWeight = 'bold';
    ax.FontName = 'Cambria';
    ax.FontSize = font_size;

end

