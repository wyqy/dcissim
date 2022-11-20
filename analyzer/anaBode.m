function handle_bode = anaBode(ss_original, ss_identified)
%ANABODE 绘制Bode图
    
    % 参数
    orilegend = 'Original';
    indlegend = 'Identified';
    font_size = 10;
    y_size = size(ss_original.C, 1);
    u_size = size(ss_original.B, 2);

    % 图窗
    figure;
    % f = figure; f.Units = 'centimeters'; f.Position = [0 0 14 10];
    
    % 变量预定义
    state_cell_row = cell(u_size, 1);  % 行表示输入
    state_cell_col = cell(y_size, 1);  % 列表示输出
    for iter_init = 1:u_size, state_cell_row{iter_init} = 'off'; end
    for iter_init = 1:y_size, state_cell_col{iter_init} = 'off'; end
    
    for iter_ax = 1:y_size*u_size
        % 子图位置
        ax = subplot(u_size, y_size, iter_ax);
        % 参数计算
        location_row = ceil(iter_ax/y_size);  % 同时是输入标号
        locaiton_col = mod(iter_ax-1, y_size)+1;  % 同时是输出标号

        % Bode图配置
        plotoptions = bodeoptions('cstprefs');  % ctrlpref
        plotoptions.PhaseWrapping = 'on';
        plotoptions.IOGrouping = 'all';
        plotoptions.XLabel.String = 'Frequency'; plotoptions.XLabel.FontSize = font_size; plotoptions.XLabel.FontWeight = 'bold'; plotoptions.XLabel.Interpreter = 'latex';
        plotoptions.YLabel.String = {'Amplitude', 'Phase'}; plotoptions.YLabel.FontSize = font_size; plotoptions.YLabel.FontWeight = 'bold'; plotoptions.YLabel.Interpreter = 'latex';
        plotoptions.Title.String = ['u' int2str(location_row) ' - y' int2str(locaiton_col)]; plotoptions.Title.FontSize = font_size; plotoptions.Title.FontWeight = 'bold'; plotoptions.Title.Interpreter = 'latex';
        temp_state = state_cell_row; temp_state{location_row} = 'on'; plotoptions.InputVisible = temp_state;
        temp_state = state_cell_col; temp_state{locaiton_col} = 'on'; plotoptions.OutputVisible = temp_state;
        % 绘制Bode图
        % handle_bode = bodeplot(ax, ss_original, 'b-', ss_identified, 'r-', {0.1, 10}, plotoptions);
        handle_bode = bodeplot(ax, ss_original, 'b-', ss_identified, 'r-', plotoptions);
        % 在第一个图上显示图例
        if iter_ax == 1, bodeLabelCaller(orilegend, indlegend, font_size); end
    end
    
    % 保存图片
    % filepath = 'fig';
    % saveas(gcf, [filepath 'steering_lateral_bode_plot.eps'], 'epsc')
end

function bodeLabelCaller(text1, text2, font_size)
    lgd = legend(text1, text2, 'FontName', 'Times New Roman', 'FontSize', font_size, 'FontWeight', 'bold');
    lgd.Location = 'northwest'; lgd.Box = 'off';
end


