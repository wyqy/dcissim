function [norm2, norminf] = anaNorm(ss_model_or_matrix)
%ANANORM 计算离散系统的H_2与H_\infty范数, 或者系统矩阵的F与2范数

    if isnumeric(ss_model_or_matrix)
        norm2 = norm(ss_model_or_matrix, 'fro');
        norminf = norm(ss_model_or_matrix, 2);
    else
        norm2 = norm(ss_model_or_matrix, 2);
        norminf = hinfnorm(ss_model_or_matrix);
    end
    
    
end

