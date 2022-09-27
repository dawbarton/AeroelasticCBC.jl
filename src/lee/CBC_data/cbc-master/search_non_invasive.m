
% initial_freq=freq_actual;
initial_freq=2.4646;
% initial_amp=amplitude/5*3;
% mean_pos=3.3647;

initial_amp=0.06;
rtc.par.x1_target_coeffs(idx_sin)=initial_amp;
% define the target A0 + A1 sin(2pi f t)
rtc.par.x1_target_coeffs(9)=mean_pos; % A0 %for the constant term in the fourier expansion
rtc.par.forcing_freq= initial_freq;

rtc.par.x_Kp=15;
% Effect of Kd:
rtc.par.x_Kd =-0.05;

error1=100;
rtc.par.x1_control =1;
search_index=1;

while sqrt(abs(error1))>error_tol

    picard=1;
    
    exp = cbc_create_experiment(rtc);
    [success,freq_actual,amplitude,mean_pos] = cbc_wait_for_convergence2(exp);

    error1=amplitude^2-rtc.par.x1_target_coeffs(idx_sin)^2;

    search_point(search_index,:)=[rtc.par.x1_target_coeffs(idx_sin);error1];
    
    if ~success
        pt = cbc_get_data_point(exp);
        failed_gain=pt;
        error('%s%s','Hopf Bifurcation: Failed to converge to steady-state initially', datestr(now, 13));
    end
    
    % Picard in the mean A0. At this point the frequency is not yet updated.
    % You could perform an iteration also in the frequency at this stage if
    % needed. The aim is just to see if the control can reduce the error in amplitude once
    % the target is fixed. Moreover, given the phase projection adopted in the
    % PP-CBC, the frequency should be about the same of the target.
    
    while abs( mean_pos-rtc.par.x1_target_coeffs(9))/abs(rtc.par.x1_target_coeffs(9))>PicardTol && abs(freq_actual-rtc.par.forcing_freq)>PicardTol_Freq
        picard
        freq_error=abs(freq_actual-rtc.par.forcing_freq);
        
        rtc.par.x1_target_coeffs(9)=mean_pos; % needed to have then the sine wave shift and the control around zero
        rtc.par.forcing_freq= freq_actual;
        
        exp = cbc_create_experiment(rtc);
        [success,freq_actual,amplitude,mean_pos] = cbc_wait_for_convergence2(exp);
        
        if ~success
            pt = cbc_get_data_point(exp);
            failed_gain=pt;
            error('%s%s','Hopf Bifurcation: Failed to converge to steady-state initially', datestr(now, 13));
        end
        
        picard=picard+1;
    end
    error1=amplitude^2-rtc.par.x1_target_coeffs(idx_sin)^2;
    error2=sqrt(abs(error1))
    search_index=search_index+1;
    search_point(search_index,:)=[rtc.par.x1_target_coeffs(idx_sin);error1];
    sign_search=search_point(search_index,2)*search_point(search_index-1,2)

    if sign_search>0        
        if error1>0
            rtc.par.x1_target_coeffs(idx_sin)=rtc.par.x1_target_coeffs(idx_sin)-0.005;
        else
            rtc.par.x1_target_coeffs(idx_sin)=rtc.par.x1_target_coeffs(idx_sin)+0.005;
        end
    
    else
        rtc.par.x1_target_coeffs(idx_sin)=(search_point(search_index,1)+search_point(search_index-1,1))/2;
    end
end

n = 100000; % sampling points
skip = 1; %
time = (1:n)/(double(rtc.par.sample_freq)/skip);
rtc.set_stream(0, {'x1', 'x1_target','Fshaker','aksim_angle','mean_h'}, n,0);
data = rtc.run_stream(0);
% save('CBC_unstalbe_v15_6_2')

plot(data(1,:))
hold on
plot(data(2,:))