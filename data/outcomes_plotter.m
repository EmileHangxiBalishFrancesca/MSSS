close all

N_p_vector = [36];   % number of people
N_t_vector = [6:2:14];  % number of table % should be ven if it is rectangular % 
tableShape_vector = [1:2]; % 1 for circle 2 for rectangular
dist_f_vector = [1]; %distance food-food


trial=0;
for i=1:length(N_p_vector)
    %for j=1:length(N_t_vector)
        for k=1:length(tableShape_vector)
            for l=1:length(dist_f_vector)
                trial = trial+1;
                hold on
                plot(N_t_vector,cost_time_tot_matrix(i,:,k,l))
                hold off
                string1 = strcat('N_p=',num2str(N_p_vector(i),'%.3d'));
                %string2 = strcat('   N_t=',num2str(N_t_vector(j),'%.3d'));
                string3 = strcat('   table_shape=',num2str(tableShape_vector(k),'%.3d'));
                string4 = strcat('   dist_f=',num2str(dist_f_vector(l),'%.3d'));

                legend_sentence(trial,:) = strcat(string1,string3,string4);
                
            end
        end
    end
%end
legend(legend_sentence)