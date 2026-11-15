function [x11, x12, x21, x22] = banlance_fast(A1, A2, b1, b2, K11, K12, K21, K22, x_old)
    % A1: New_Ag, A2: New_Ag_decay
    % b1: Ig_tot, b2: Ig2_tot
    % K11/K12: Ig对主/亚表位亲和力; K21/K22: Ig2对主/亚表位亲和力
    % x_old: [x11; x12; x21; x22] 上一时刻解，用于热启动

    if nargin < 9 || isempty(x_old)
        X = [1e-12; 1e-12; 1e-12; 1e-12];
    else
        X = x_old;
    end

    tol = 1e-13;
    max_iter = 20;
    invK = [1/K11; 1/K12; 1/K21; 1/K22];

    for iter = 1:max_iter
        % 提取当前分量
        x11=X(1); x12=X(2); x21=X(3); x22=X(4);

        % 计算游离量 (Free amounts)
        phi1 = A1 - x11 - x21; % Free New_Ag
        phi2 = A2 - x12 - x22; % Free New_Ag_decay
        psi1 = b1 - x11 - x12; % Free Ig
        psi2 = b2 - x21 - x22; % Free Ig2

        % 方程残差 F
        F = [phi1 * psi1 - x11 * invK(1);
             phi2 * psi1 - x12 * invK(2);
             phi1 * psi2 - x21 * invK(3);
             phi2 * psi2 - x22 * invK(4)];

        if norm(F, inf) < tol, break; end

        % 构造雅可比矩阵 J (4x4)
        % 行1对应F1对x11,x12,x21,x22的导数
        J = [ -(psi1+phi1+invK(1)), -phi1,               -psi1,                0;
              -phi2,               -(psi1+phi2+invK(2)),  0,                  -psi1;
              -psi2,                0,               -(psi2+phi1+invK(3)), -phi1;
               0,                  -psi2,               -phi2,               -(psi2+phi2+invK(4))];

        % 更新步长
        dX = J \ (-F);
        X = X + dX;

        % 物理约束边界
        X = max(1e-18, X); 
        X(1) = min(X(1), min(A1, b1)*0.999);
        X(2) = min(X(2), min(A2, b1)*0.999);
        X(3) = min(X(3), min(A1, b2)*0.999);
        X(4) = min(X(4), min(A2, b2)*0.999);
    end
    x11=X(1); x12=X(2); x21=X(3); x22=X(4);
end