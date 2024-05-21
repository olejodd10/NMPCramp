function TestMexLtiMpc(input_dir, output_dir, N, simulation_timesteps)
    % Input data
    % Matrices are transposed to be compatible with C
    Q = csvread(input_dir + "/Q.csv")';
    S = csvread(input_dir + "/S.csv")';
    R = csvread(input_dir + "/R.csv")';
    fx = csvread(input_dir + "/fx.csv")';
    fu = csvread(input_dir + "/fu.csv")';
    A = csvread(input_dir + "/A.csv")';
    B = csvread(input_dir + "/B.csv")';
    C = csvread(input_dir + "/C.csv")';
    y_min = csvread(input_dir + "/y_min.csv")';
    y_max = csvread(input_dir + "/y_max.csv")';
    Lt = csvread(input_dir + "/_Lt.csv")';
    lt = csvread(input_dir + "/lt.csv")';
    u_min = csvread(input_dir + "/u_min.csv")';
    u_max = csvread(input_dir + "/u_max.csv")';
    x0 = csvread(input_dir + "/x0.csv")';

    n_x = size(Q, 1);
    n_u = size(R, 1);
    n_y = size(C, 2);
    n_t = size(Lt, 2);

    % Pretend system is time varying
    A = repmat(A, 1, 1, N);
    B = repmat(B, 1, 1, N);
    d = zeros(n_x, N, 'double');
    fx = repmat(fx', 1, N);
    fu = repmat(fu', 1, N);

    % State and input allocations, and initial state
    xout = zeros(simulation_timesteps+1, n_x, 'double');
    xout(1,:) = x0;
    uout = zeros(simulation_timesteps, n_u, 'double');

    % Simulate
    simulation_time_s = 0.0;
    for i = 1:simulation_timesteps
        tic;
        % Solve QP
        xout(i,:); % Some weird precision bug makes this necessary
        [x,u] = MexLtvMpc(n_x, n_u, n_y, n_t, N, Q, S, R, fx, fu, A, B, d, C, y_min, y_max, Lt, lt, u_min, u_max, xout(i,:));
        simulation_time_s = simulation_time_s + toc;

        % Simulate using linearized discrete model
        u0 = u(:,1);
        xout(i+1,:) = (A(:,:,1)'*xout(i,:)' + B(:,:,1)'*u0)';
        uout(i,:) = u(:,1)';
    end

    % Timer output
    fprintf("%d timesteps with horizon %d finished in %f ms\n", simulation_timesteps, N, simulation_time_s*1000);

    % Save output
    writematrix(xout, output_dir + "/xout.csv");
    writematrix(uout, output_dir + "/uout.csv");

    clear MexLtvMpc;
end
