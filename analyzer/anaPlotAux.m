function anaPlotAux(obj, font_size, line_width)
%ANAPLOTAUX 修改axis的格式

    % 常规美化
    if isgraphics(obj, 'Axes')
        obj.LineWidth = line_width;
        obj.Box = 'on';
        obj.TickDir = 'in';
        obj.XMinorTick = 'on';
        obj.YMinorTick = 'on';
        obj.XLimitMethod = 'tight';
        obj.XGrid = 'on';
        obj.YGrid = 'on';
        obj.GridLineStyle = '--';
        obj.XColor = [.3,.3,.3];
        obj.YColor = [.3,.3,.3];
        obj.FontWeight = 'bold';
        obj.FontName = 'Arial';
        obj.FontSize = font_size;
    elseif isgraphics(obj, 'Line')
        obj.LineWidth = line_width;
    end

end

