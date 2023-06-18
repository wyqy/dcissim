function [err_h2, err_hinf] = errSystem(sys_ori, sys_ind)
%ERRSYSTEM 比较系统的H2和H-\infty范数

    [num_h2, num_hinf] = anaNorm(sys_ind - sys_ori);
    [den_h2, den_hinf] = anaNorm(sys_ori);
    err_h2 = abs(num_h2)/den_h2;
    err_hinf = abs(num_hinf)/den_hinf;
end


function [norm2, norminf] = anaNorm(ss_model_or_matrix)
    % 计算离散系统的H_2与H_\infty范数

    norm2 = norm(ss_model_or_matrix, 2);
    norminf = hinfnorm(ss_model_or_matrix);
    
end
