close all
clear legend_sentence
%THE END OF THE LOOPS
N_p_vector = [72];   % number of people
N_t_vector = [8];  % number of table % should be ven if it is rectangular %
tableShape_vector = [1:2]; % 1 for circle 2 for rectangular
dist_f_vector = [0.01:0.3:3.0]; %distance food-food
number_statistical_attemps = 20;


trial=0;
for i=1:length(N_p_vector)
    for j=1:length(N_t_vector)
        for k=1:length(tableShape_vector)
            %for l=1:length(dist_f_vector)
                trial = trial+1;
                figure(1)
                hold on
                errorbar(dist_f_vector,reshape(cost_time_tot_matrix(i,j,k,:),1,[]),reshape(variance_time(i,j,k,:),1,[]))
                
                hold off
                figure(2)
                hold on
                errorbar(dist_f_vector,reshape(cost_f_tot_matrix(i,j,k,:),1,[]),reshape(variance_f(i,j,k,:),1,[]))
                
                hold off
                figure(3)
                hold on
                errorbar(dist_f_vector,reshape(cost_v_tot_matrix(i,j,k,:),1,[]),reshape(variance_v(i,j,k,:),1,[]))
                
                hold off
                
                string1 = strcat('N_p=',num2str(N_p_vector(i),'%.3d'));
                string2 = strcat('   N_t=',num2str(N_t_vector(j),'%.3d'));
                string3 = strcat('   table_shape=',num2str(tableShape_vector(k),'%.3d'));
                %string4 = strcat('   dist_f=',num2str(dist_f_vector(l),'%.3d'));

                legend_sentence(trial,:) = strcat(string1,string2,string3);
                
            end
        end
    end
%end

clear legend_sentence
legend_sentence = ['Circular table disposition   ';'Rectangular table disposition'];

figure(1)
legend(legend_sentence,'Location','northwest')
title('Cost related to the time')
xlabel('Distance between food locations')
ylabel('Velocity cost')
figure(2)
legend(legend_sentence,'Location','northwest')
title('Cost related to the forces')
xlabel('Distance between food locations')
ylabel('Velocity cost')
figure(3)
legend(legend_sentence,'Location','northwest')
title('Cost related to the velocity')
xlabel('Distance between food locations')
ylabel('Velocity cost')
