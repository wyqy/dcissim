function [y_isim, u_isim] = idenISIM(yn, un, vn, isim_lsq_type)
% IDENISIM 辨识U(v)和Y(v)矩阵

    % recursive变量数值存储在这里
    persistent is_inited y_size u_size v_size mat_lambda mat_p mat_k vec_r

    % 不同辨识方法 (选择一种即可)
    switch isim_lsq_type
        case 'ordinary' % 离线: 直接做QR分解和COD分解的OLS
            % 参数计算
            y_size = size(yn, 1); u_size = size(un, 1); v_size = size(vn, 1);
            signal_sample = size(yn, 2);
            yn = yn.'; un = un.'; vn = vn.';
            yun = [yn un];
            % 最小二乘
            % yu_isim = lsqminnorm(vn, yun);
            [yu_isim, flag] = lsqr(kron(eye(y_size+u_size), vn), reshape(yun, [(y_size+u_size)*signal_sample 1]), 1e-8);
            if flag ~= 0, yu_isim = zeros(v_size, (y_size+u_size));
            else, yu_isim = reshape(yu_isim, [v_size (y_size+u_size)]); end
            % 返回值
            yu_isim = yu_isim.';
            y_isim = yu_isim(1:y_size, :);
            u_isim = yu_isim(y_size+1:end, :);
            
        case 'recursive' % 在线: 逐个样本点做RLS
            % 首次使用, 初始化
            if isempty(is_inited)
                % 参数计算
                y_size = size(yn, 1); u_size = size(un, 1); v_size = size(vn, 1);
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
                vec_yu = [yn; un]; vec_v = vn;
                mat_phi = kron(vec_v, eye(y_size+u_size));
                
                dmat_a = decomposition(mat_lambda + (mat_phi.') * mat_p * mat_phi);
                mat_k = mat_p * mat_phi / dmat_a;
                mat_p = mat_p - mat_p * mat_phi / dmat_a * (mat_phi.') * mat_p;
                vec_r = vec_r + mat_k * (vec_yu - (mat_phi.') * vec_r);

                % 返回数值
                yu_isim = reshape(vec_r, [(y_size+u_size) v_size]);
                y_isim = yu_isim(1:y_size, :);
                u_isim = yu_isim(y_size+1:end, :);
            end
        otherwise, y_isim = 0; u_isim = 0;
    end

end
