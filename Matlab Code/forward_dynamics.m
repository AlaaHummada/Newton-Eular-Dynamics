function [theta, theta_dot] = forward_dynamics(theta0, theta_dot0, tau, F_tip, M, A, G, g, tf, N)
    % Inputs:
    %   theta0: Initial joint positions (n x 1 vector)
    %   theta_dot0: Initial joint velocities (n x 1 vector)
    %   tau: Joint torques (n x 1 vector or function handle tau(t))
    %   F_tip: External wrench at end-effector (6 x 1 vector or function handle F_tip(t))
    %   M: Cell array of M_{i,i-1} configuration matrices
    %   A: Cell array of screw axes A_i for each joint
    %   G: Cell array of spatial inertia matrices G_i for each link
    %   g: Gravity vector in base frame (3 x 1 vector)
    %   tf: Final time
    %   N: Number of integration steps
    %
    % Outputs:
    %   theta: Joint position trajectory (n x N+1 matrix)
    %   theta_dot: Joint velocity trajectory (n x N+1 matrix)

    n = length(theta0);
    dt = tf / N;
    
    % Initialize trajectories
    theta = zeros(n, N+1);
    theta_dot = zeros(n, N+1);
    theta(:,1) = theta0;
    theta_dot(:,1) = theta_dot0;

    % Euler Integration
    for k = 1:N
        t = (k-1) * dt;
        
        % Evaluate tau and F_tip at time t
        if isa(tau, 'function_handle')
            tau_k = tau(t);
        else
            tau_k = tau;
        end
        if isa(F_tip, 'function_handle')
            F_tip_k = F_tip(t);
        else
            F_tip_k = F_tip;
        end
        
        % Compute joint acceleration
        theta_ddot = compute_forward_dynamics(theta(:,k), theta_dot(:,k), tau_k, F_tip_k, M, A, G, g);
        
        % Update states using Euler integration
        theta(:,k+1) = theta(:,k) + theta_dot(:,k) * dt;
        theta_dot(:,k+1) = theta_dot(:,k) + theta_ddot * dt;
    end
end

function theta_ddot = compute_forward_dynamics(theta, theta_dot, tau, F_tip, M, A, G, g)
    % Helper function to compute joint accelerations
    n = length(theta);
    
    % Compute h(theta, theta_dot)
    h = inverse_dynamics(theta, theta_dot, zeros(n,1), zeros(6,1), M, A, G, g);
    
    % Compute mass matrix M(theta)
    M_theta = zeros(n, n);
    for i = 1:n
        theta_ddot_i = zeros(n,1);
        theta_ddot_i(i) = 1;
        M_theta(:,i) = inverse_dynamics(theta, zeros(n,1), theta_ddot_i, zeros(6,1), M, A, G, zeros(3,1));
    end
    
    % Compute Jacobian J(theta) (assumed to be provided or computed separately)
    J = compute_jacobian(theta, M, A); % Placeholder for Jacobian computation
    
    % Solve for theta_ddot
    theta_ddot = M_theta \ (tau - h - J' * F_tip);
end

function J = compute_jacobian(theta, M, A)
    n = length(theta);
    J_s = zeros(6, n);
    M_0_i = cell(n, 1);
    T = eye(4);
    for i = 1:n
        T = T * inv(M{i});
        M_0_i{i} = T;
    end
    S = cell(n, 1);
    for i = 1:n
        Ad_M_0_i = adjoint(M_0_i{i});
        S{i} = Ad_M_0_i * A{i};
    end
    T = eye(4);
    for i = 1:n
        Ad_T = adjoint(T);
        J_s(:,i) = Ad_T * S{i};
        T = T * expm(skew(S{i}) * theta(i));
    end
    Ad_T_inv = adjoint(inv(T));
    J = Ad_T_inv * J_s;
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
    R = T(1:3, 1:3);
    p = T(1:3, 4);
    Ad = [R, zeros(3); skew(p) * R, R];
end

function ad = ad_matrix(V)
    omega = V(1:3);
    v = V(4:6);
    ad = [skew(omega), zeros(3); skew(v), skew(omega)];
end