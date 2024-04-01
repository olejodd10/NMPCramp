function TestMexSimulateMmc(output_dir, N, simulation_timesteps)
    % System parameters
    R = 10.0e-3;
    L = 1.5e-3;
    C = 20.0e-3;
    
    Rc = 0.0;
    Lc = 0.0;
    
    N_SM = 18.0;
    
    FREQ = 50.0;
    PHASE_0 = 0.0;
    
    Vf_amp = 100.0;
    Vf_phase = 0.0; % Phase relative to PHASE_0
    Vdc = 300.0;
    
    % Costs and references
    q1 = 1.0;
    q2 = 0.3;
    
    Iv_ref_amp = 50.0;
    Iv_phase = pi; % Phase relative to PHASE_0
    
    P = 7.5e3;
    Idc_ref = P/Vdc;
    Icir_ref = Idc_ref/3.0;
    
    % Model and discretization
    N_X = 4;
    N_U = 2;
    
    Ts = 70.0e-6;
    
    % Constraints
    Iv_0_min = -75.0;
    Icir_0_min = -50.0;
    Vsigma_u_min = 0.0;
    Vsigma_l_min = 0.0;
    
    X_MIN = [Iv_0_min, Icir_0_min, Vsigma_u_min, Vsigma_l_min];
    
    Iv_0_max = 75.0;
    Icir_0_max = 50.0;
    Vsigma_u_max = 400.0;
    Vsigma_l_max = 400.0;
    
    X_MAX = [Iv_0_max, Icir_0_max, Vsigma_u_max, Vsigma_l_max];
    
    INSERTION_INDEX_DEVIATION_ALLOWANCE = 2.0;
    
    u1_min = 0.0;
    u2_min = 0.0;
    
    U_MIN = [u1_min, u2_min];
    
    u1_max = N_SM;
    u2_max = N_SM;
    
    U_MAX = [u1_max, u2_max];
    
    % Initial conditions
    Iv_0 = Iv_ref_amp*sin(PHASE_0 + Iv_phase);
    Icir_0 = Icir_ref;
    Vsigma_u_0 = Vdc;
    Vsigma_l_0 = Vdc;

    % State and input allocations, and initial state
    xout = zeros(simulation_timesteps+1, N_X, 'double');
    xout(1,:) = [Iv_0, Icir_0, Vsigma_u_0, Vsigma_l_0];
    uout = zeros(simulation_timesteps, N_U, 'double');

    % Allocate and initialize trajectory
    % This will be used by C code, therefore the dimensions are transposed
    x = zeros(N_X, N, 'double');
    u = zeros(N_U, N, 'double');
    for i = 1:N
        x(1,i) = Iv_ref_amp*sin(PHASE_0 + Iv_phase + 2*pi*FREQ*Ts*(i-1));
        x(2,i) = Icir_0;
        x(3,i) = Vsigma_u_0;
        x(4,i) = Vsigma_l_0;
        u(1,i) = (N_SM/2)*sin(PHASE_0 - Vf_phase + 2*pi*FREQ*Ts*(i-1)) + N_SM/2;
        u(2,i) = (N_SM/2)*sin(PHASE_0 + Vf_phase + 2*pi*FREQ*Ts*(i-1)) + N_SM/2;
    end

    % Allocate model matrices
    A = zeros(N_X, N_X, N, 'double');
    B = zeros(N_U, N_X, N, 'double');
    d = zeros(N_X, N, 'double');

    % Simulate
    for i = 1:simulation_timesteps
        % Predict trajectory
        x(:,1) = xout(i,:)'; % We prefer simulation/measurement value to MPC value
        % Shift input trajectory
        u = [u(:,2:end), u(:,end)];

        % Extrapolate references and disturbances
        vf = Vf_amp*sin(PHASE_0 + Vf_phase + 2*pi*FREQ*Ts*(i-1 + (0:(N-1))));
        Iv_ref = Iv_ref_amp*sin(PHASE_0 + Iv_phase + 2*pi*FREQ*Ts*(i-1 + (0:(N-1))));

        % Get linearized discrete model
        [A,B,d] = MexMmcModel(R, Rc, L, Lc, C, Ts, N_SM, N, x, u, vf, Vdc, A, B, d);

        % Solve QP
        [x,u] = MexSdqpLmpcMmc(N_X, N_U, N, q1, q2, X_MIN, X_MAX, N_SM, INSERTION_INDEX_DEVIATION_ALLOWANCE, U_MIN, U_MAX, Iv_ref, Icir_ref, A, B, d, xout(i,:), x, u);

        % Simulate using linearized discrete model
        A0 = squeeze(A(:,:,1))';
        B0 = squeeze(B(:,:,1))';
        d0 = squeeze(d(:,1));
        x0 = xout(i,:)';
        u0 = u(:,1);
        xout(i+1,:) = (A0*x0 + B0*u0 + d0)';
        uout(i,:) = u(:,1)';
    end

    % Save output
    writematrix(xout, output_dir + "/xoutN" + N + ".csv");
    writematrix(uout, output_dir + "/uoutN" + N + ".csv");

end
