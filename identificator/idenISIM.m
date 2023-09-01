function [yv_est, uv_est, yw_est, uw_est] = idenISIM(yk, uk, mat_v, T, isim_lsq_type)
% IDENISIM 辨识U(v)和Y(v)矩阵

    % recursive变量数值存储在这里
    persistent is_inited y_size u_size v_size mat_w mat_r k t

    % 不同辨识方法 (选择一种即可)
    switch isim_lsq_type
        case 'ols' % 离线
            % 参数计算
            y_size = size(yk, 1); u_size = size(uk, 1); v_size = size(mat_v, 1);
            N = size(yk, 2); P = floor(N/T);

            % 矩阵转置
            yk = yk.';      uk = uk.';
            yuk = [yk uk];
            yuk = yuk(1:P*T, :);
            % 按周期分块相加求和
            yuk = reshape(yuk, [T, P, y_size+u_size]);
            yuk = squeeze(sum(yuk, 2));
            % 乘法
            yuv_est = (2/(P*T)) .* mat_v * yuk;
            
            % 返回值转置
            yuv_est = yuv_est.';
            yv_est = yuv_est(1:y_size, :);      uv_est = yuv_est(y_size+1:end, :);
            yw_est = 0; uw_est = 0;

        case 'fft' % FFT
            % 参数计算
            y_size = size(yk, 1);   u_size = size(uk, 1);   v_size = size(mat_v, 1);
            N = size(yk, 2);        P = floor(N/T);
            
            % 输入值准备
            yk = yk(:, 1:P*T);    uk = uk(:, 1:P*T);
            % 返回值准备
            yw_est = zeros(y_size, T);          uw_est = zeros(u_size, T);
            yv_est = zeros(y_size, v_size);     uv_est = zeros(u_size, v_size);

            % 按周期分块FFT
            for iter_k = 0:P-1
                yw_est = yw_est + fft(yk(:, iter_k*T+1:iter_k*T+T), T, 2);
                uw_est = uw_est + fft(uk(:, iter_k*T+1:iter_k*T+T), T, 2);
            end
            yw_est = yw_est ./ N; uw_est = uw_est ./ N;

            % 变换
            temp_list = 2:2:v_size; temp_size = length(temp_list);
            yv_est(:, temp_list) = -2 .* imag(yw_est(:, 2:temp_size+1));    uv_est(:, temp_list) = -2 .* imag(uw_est(:, 2:temp_size+1));
            temp_list = 3:2:v_size; temp_size = length(temp_list);
            yv_est(:, temp_list) = 2 .* real(yw_est(:, 2:temp_size+1));     uv_est(:, temp_list) = 2 .* real(uw_est(:, 2:temp_size+1));
            yv_est(:, 1) = sqrt(2) .* real(yw_est(:, 1));                   uv_est(:, 1) = sqrt(2) .* real(uw_est(:, 1));
            
        case 'rls' % 在线
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
                yv_est = yuv_est(1:y_size, :); uv_est = yuv_est(y_size+1:end, :);
                yw_est = 0; uw_est = 0;
            end

        otherwise, yv_est = 0; uv_est = 0; yw_est = 0; uw_est = 0;
    end


end
