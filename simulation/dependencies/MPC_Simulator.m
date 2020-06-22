%Signal Based MPC Function for Fixed-Point Analysis


function [y_t, u_t,quadprog_runtimes, pd_runtimes] = MPC_Simulator(ref_1, ref_2, T,H,q_,E,F1,z_min,z_max,n,m,N_h,A_d,B_d,C_d,Aob,Bob,Cob,Ld) %#codegen


    % -------- System Initial Conditions ---------

    x_d_obs = zeros(n+m,T); 
    x_plant = zeros(n,T); 
    u_t = zeros(m,T);
    z = zeros( (n+m)*N_h, T);
    y_t = zeros(m,T); 
    y_obs = zeros(m,T); 
    
    quadprog_runtimes = [];
    pd_runtimes = [];
    
    for i = 1:T-1

        %Simulate Plant Response
        x_plant(:,i+1) = A_d*x_plant(:,i) + B_d*u_t(:,i); 
        y_t(:,i) = C_d * x_plant(:,i);


        %Observer 
        x_d_obs(:,i+1) = Aob * x_d_obs(:,i) + Bob*u_t(:,i) - Ld * ( Cob*x_d_obs(:,i) - y_t(:,i) );


        %Scale Reference Voltage to height 
        curr_ref_t0 = ref_1(:,i) * 7.15; %V to Cm scaling factor
        curr_ref_t1 = ref_2(:,i) * 7.15; 

        curr_ref = [curr_ref_t0; curr_ref_t1];

        %Compute f and e Matrix 
        f = q_ * (curr_ref - x_d_obs(n+1:end,i)); 
        e = F1 * x_d_obs(1:n,i);  
        
        % Benchmark against quadprog execution time
        tic;
        U= quadprog(H,f,[],[],[],[],z_min,z_max, []);
        quadprog_time = toc; 
        guadprog_runtimes = [quadprog_runtimes quadprog_time]; 
        
        tic; 
        z(:,i) = PrimalDual(H,-f,E,e,z_min,z_max);    
        primaldual_time = toc; 
        pd_runtimes = [pd_runtimes primaldual_time]; 

        %Implement DAC Based Saturation
        for j = 1:2 
            if z(j,i) < -3.3 
                u_t(j,i+1) = z(j,i); 
            elseif z(j,i) > 3.3 
                u_t(j,i+1) = z(j,i); 
            else     
                u_t(j,i+1) = z(j,i); 

            end 
        end 

    end 

end 