close all
figure(1)
hold on
%{
for i=1:20
    plot(1:375,cost_f_mean(:,i))
end
%}


cost_f_tot_matrix_norm = cost_f_tot_matrix / max(cost_f_tot_matrix(:));
cost_v_tot_matrix_norm = cost_v_tot_matrix / max(cost_v_tot_matrix(:));
cost_time_tot_matrix_norm = cost_time_tot_matrix / max(cost_time_tot_matrix(:));


p1 = 13;
p2 = 7;

cost_f_tot_matrix_norm = reshape(cost_f_tot_matrix_norm(:,1,2,:),p1,p2);
cost_v_tot_matrix_norm = reshape(cost_v_tot_matrix_norm(:,1,2,:),p1,p2);
cost_time_tot_matrix_norm = reshape(cost_time_tot_matrix_norm(:,1,2,:),p1,p2);

[X,Y] = meshgrid(10:2:34,0.2:0.2:1.4);

scatter3(X(:),Y(:),cost_f_tot_matrix_norm(:),'ro')
scatter3(X(:),Y(:),cost_v_tot_matrix_norm(:),'bo')
scatter3(X(:),Y(:),cost_time_tot_matrix_norm(:),'go')



hold off
