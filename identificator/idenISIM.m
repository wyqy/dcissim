function [yv_est, uv_est] = idenISIM(yk, uk, vk, isim_lsq_type)
% IDENISIM 辨识U(v)和Y(v)矩阵

    % recursive变量数值存储在这里
    persistent is_inited y_size u_size v_size mat_lambda mat_p mat_k vec_r

    % 不同辨识方法 (选择一种即可)
    switch isim_lsq_type
        case 'ordinary' % 离线: 直接做QR分解和COD分解的OLS
            % 参数计算
            y_size = size(yk, 1);
            u_size = size(uk, 1);
            v_size = size(vk, 1);
            sample_size = size(yk, 2);
            % 准备
            yk = yk.';
            uk = uk.';
            vk = vk.';
            yuk = [yk uk];

            % 最小二乘
            if v_size > sample_size/2  % 频率较多
                [yuv_est, flag] = lsqr(kron(eye(y_size+u_size), vk), reshape(yuk, [(y_size+u_size)*sample_size 1]), 1e-6, 100);
                if flag ~= 0, yuv_est = zeros(v_size, (y_size+u_size));
                else, yuv_est = reshape(yuv_est, [v_size (y_size+u_size)]); end
            else  % 频率较少
                yuv_est = lsqminnorm(vk, yuk);
            end
            % 返回值
            yuv_est = yuv_est.';
            yv_est = yuv_est(1:y_size, :);
            uv_est = yuv_est(y_size+1:end, :);
            
        case 'recursive' % 在线: 逐个样本点做RLS
            % 首次使用, 初始化
            if isempty(is_inited)
                % 参数计算
                y_size = size(yk, 1); u_size = size(uk, 1); v_size = size(vk, 1);
                % 矩阵生成
                seed = rng().Seed;
                rs = RandStream('dsfmt19937', 'Seed', seed);
                mat_lambda = eye(y_size+u_size);
                mat_p = 10e5 * eye((y_size+u_size)*v_size) + randi(rs, 10e4);
                vec_r = zeros((y_size+u_size)*v_size, 1);
                % 设置标记
                is_inited = 1;
            % 已有值, 迭代一步
            else
                % 开始迭代
                vec_yu = [yk; uk]; vec_v = vk;
                mat_phi = kron(vec_v, eye(y_size+u_size));
                
                dmat_a = decomposition(mat_lambda + (mat_phi.') * mat_p * mat_phi);
                mat_k = mat_p * mat_phi / dmat_a;
                mat_p = mat_p - mat_p * mat_phi / dmat_a * (mat_phi.') * mat_p;
                vec_r = vec_r + mat_k * (vec_yu - (mat_phi.') * vec_r);

                % 返回数值
                yuv_est = reshape(vec_r, [(y_size+u_size) v_size]);
                yv_est = yuv_est(1:y_size, :);
                uv_est = yuv_est(y_size+1:end, :);
            end
        otherwise, yv_est = 0; uv_est = 0;
    end

end
