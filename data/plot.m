%load('CHANGING--N_p--tableShape--dist_f--TIME_15-Dec-2017 12:32:08.mat')
%time
%distance 0.5
subplot(2,2,1)
plot(N_p_vector,cost_time_tot_matrix(:,:,1,1),'b--o') 
hold on
plot(N_p_vector,cost_time_tot_matrix(:,:,2,1),'b--*')

%distance 1.5
plot(N_p_vector,cost_time_tot_matrix(:,:,1,2),'c--o') 
hold on
plot(N_p_vector,cost_time_tot_matrix(:,:,2,2),'c--*')

%distance 2.5
plot(N_p_vector,cost_time_tot_matrix(:,:,1,3),'g--o') 
hold on
plot(N_p_vector,cost_time_tot_matrix(:,:,2,3),'g--*')

%distance 3.5
plot(N_p_vector,cost_time_tot_matrix(:,:,1,4),'r--o') 
hold on
plot(N_p_vector,cost_time_tot_matrix(:,:,2,4),'r--*')
legend('circle 0.5','rectangle 0.5','circle 1.5','rectangle 1.5','circle 2.5','rectangle 2.5', 'circle 3.5','rectangle 3.5')
title('cost function TIME')
%velocity
%distance 0.5
subplot(2,2,2)
plot(N_p_vector,cost_v_tot_matrix(:,:,1,1),'b--o') 
hold on
plot(N_p_vector,cost_v_tot_matrix(:,:,2,1),'b--*')
legend('circle','rectangle')

%distance 1.5
plot(N_p_vector,cost_v_tot_matrix(:,:,1,2),'c--o') 
hold on
plot(N_p_vector,cost_v_tot_matrix(:,:,2,2),'c--*')

%distance 2.5
plot(N_p_vector,cost_v_tot_matrix(:,:,1,3),'g--o') 
hold on
plot(N_p_vector,cost_v_tot_matrix(:,:,2,3),'g--*')

%distance 3.5
plot(N_p_vector,cost_v_tot_matrix(:,:,1,4),'r--o') 
hold on
plot(N_p_vector,cost_v_tot_matrix(:,:,2,4),'r--*')
title('cost function VELOCITY')
legend('circle 0.5','rectangle 0.5','circle 1.5','rectangle 1.5','circle 2.5','rectangle 2.5', 'circle 3.5','rectangle 3.5')

%force
subplot(2,2,3)
%distance 0.5
plot(N_p_vector,cost_f_tot_matrix(:,:,1,1),'b--o') 
hold on
plot(N_p_vector,cost_f_tot_matrix(:,:,2,1),'b--*')
legend('circle','rectangle')

%distance 1.5
plot(N_p_vector,cost_f_tot_matrix(:,:,1,2),'c--o') 
hold on
plot(N_p_vector,cost_f_tot_matrix(:,:,2,2),'c--*')

%distance 2.5
plot(N_p_vector,cost_f_tot_matrix(:,:,1,3),'g--o') 
hold on
plot(N_p_vector,cost_f_tot_matrix(:,:,2,3),'g--*')

%distance 3.5
plot(N_p_vector,cost_f_tot_matrix(:,:,1,4),'r--o') 
hold on
plot(N_p_vector,cost_f_tot_matrix(:,:,2,4),'r--*')
title('cost function FORCE')
legend('circle 0.5','rectangle 0.5','circle 1.5','rectangle 1.5','circle 2.5','rectangle 2.5', 'circle 3.5','rectangle 3.5')

