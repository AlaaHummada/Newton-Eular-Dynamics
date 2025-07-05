% Include inverse_dynamics, forward_dynamics, skew, adjoint, ad_matrix functions
function tau = inverse_dynamics(theta, theta_dot, theta_ddot, F_tip, M, A, G, g)
    n = length(theta);
    nu = cell(n+1, 1);
    nu_dot = cell(n+1, 1);
    F = cell(n+1, 1);
    tau = zeros(n, 1);
    
    nu{1} = zeros(6,1);
    nu_dot{1} = [zeros(3,1); -g];
    F{n+1} = F_tip;
    
    for i = 1:n
        T_i_im1 = expm(-skew(A{i}) * theta(i)) * M{i};
        Ad_T_i_im1 = adjoint(T_i_im1);
        nu{i+1} = A{i} * theta_dot(i) + Ad_T_i_im1 * nu{i};
        ad_nu_i = ad_matrix(nu{i+1});
        nu_dot{i+1} = A{i} * theta_ddot(i) + Ad_T_i_im1 * nu_dot{i} + ad_nu_i * A{i} * theta_dot(i);
    end
    
    for i = n:-1:1
        if i == n
            T_ip1_i = M{n+1}; % For the last link, no transformation to a "next" frame
        else
            T_ip1_i = expm(-skew(A{i+1}) * theta(i+1)) * M{i+1};
        end
        Ad_T_ip1_i_T = adjoint(T_ip1_i)';
        ad_nu_i_T = ad_matrix(nu{i+1})';
        F{i} = Ad_T_ip1_i_T * F{i+1} + G{i} * nu_dot{i+1} - ad_nu_i_T * (G{i} * nu{i+1});
        tau(i) = F{i}' * A{i};
    end
end
function S = skew(v)
    % Converts a 3x1 vector to a 3x3 skew-symmetric matrix or a 6x1 vector to a 4x4 screw matrix
    if length(v) == 3
        % 3x1 vector: Return 3x3 skew-symmetric matrix
        S = [0, -v(3), v(2); v(3), 0, -v(1); -v(2), v(1), 0];
    elseif length(v) == 6
        % 6x1 vector: Return 4x4 screw matrix [S] = [[omega] v; 0 0]
        omega = v(1:3); % Angular part
        v_linear = v(4:6); % Linear part
        omega_skew = [0, -omega(3), omega(2); omega(3), 0, -omega(1); -omega(2), omega(1), 0];
        S = [omega_skew, v_linear; 0, 0, 0, 0];
    else
        error('Input vector must be 3x1 or 6x1');
    end
end

function Ad = adjoint(T)
    % Computes the 6x6 adjoint matrix for a 4x4 transformation matrix T
    R = T(1:3, 1:3);
    p = T(1:3, 4);
    Ad = [R, zeros(3); skew(p) * R, R];
end

function ad = ad_matrix(V)
    % Computes the 6x6 ad matrix for a 6x1 twist vector V
    omega = V(1:3);
    v = V(4:6);
    ad = [skew(omega), zeros(3); skew(v), skew(omega)];
end