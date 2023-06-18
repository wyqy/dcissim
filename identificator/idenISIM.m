function [yv_est, uv_est] = idenISIM(yk, uk, mat_v, T, isim_lsq_type)
% IDENISIM 辨识U(v)和Y(v)矩阵

    % recursive变量数值存储在这里
    persistent is_inited y_size u_size v_size mat_w mat_r k t

    % 不同辨识方法 (选择一种即可)
    switch isim_lsq_type
        case 'ordinary' % 离线
            % 参数计算
            y_size = size(yk, 1);
            u_size = size(uk, 1);
            v_size = size(mat_v, 1);
            N = size(yk, 2);
            
            % 输入准备
            yk = yk.';
            uk = uk.';
            yuk = [yk uk];

            % V矩阵构造
            mat_v = repmat(mat_v, [1, floor(N/T)]);

            % 直接乘法
            yuv_est = (2/(floor(N/T)*T)) .* mat_v * yuk;

            % 返回值
            yuv_est = yuv_est.';
            yv_est = yuv_est(1:y_size, :);
            uv_est = yuv_est(y_size+1:end, :);
            
        case 'recursive' % 在线
            % 首次使用, 初始化
            if isempty(is_inited)
                % 参数计算
                y_size = size(yk, 1);
                u_size = size(uk, 1);
                v_size = size(mat_v, 1);

                % 初始化
                mat_w = zeros(v_size, y_size + u_size);
                mat_r = zeros(v_size, y_size + u_size);
                k = 0; t = 0;

                % 设置标记
                is_inited = 1;
            % 已有值, 迭代一步
            else
                % 开始迭代
                k = k + 1;
                mat_w = mat_w + mat_v(:, k)*[yk.' uk.'];
                if mod(k, T) == 0
                    k = 0;
                    t = t + 1;
                    mat_w = (2/T).*mat_w;
                    mat_r = ((t-1)/t).*mat_r + (1/t).*mat_w;
                end

                % 返回数值
                yuv_est = mat_r.';
                yv_est = yuv_est(1:y_size, :);
                uv_est = yuv_est(y_size+1:end, :);
            end
        otherwise, yv_est = 0; uv_est = 0;
    end




end
