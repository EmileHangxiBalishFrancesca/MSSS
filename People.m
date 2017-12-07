close all
clear
clc
% this file is a sample program to illustract the structure of our social
% model algorithm for Apero

N_p = 30;   % number of people
tableShape = 2; % 1 for circle 2 for rectangular
N_t = 6;  % number of table % should be ven if it is rectangular
N_f = 2;   % number of foods

% Desired velocity
v0 = 0.2;
% Limit velocity
v_lim = 2*v0;

% Wall and tables constant repulsion
C_w = 0.0003;
C_t = 0.03;
% Maximum people near one table before it becomes full
max_p = 6;
% Distance at which people stay from the table centre
min_dist_table = 0.3;
% Distance at which people stay from the objective (table border of food
% point)
min_dist_obj = 0.1;
% Person-person repulsion constants
% Repulsion potential constants 
A = 2; % Taken same as the paper
B = 0.1; % same as paper

sightAngle = pi/30; % forces coming from people out of +- 60degree will be weaker
sightCoef = 0.3;

tableCapacity = 6; % capacity of 1 table can handle

%People's attributes: position, velocity, undergoing force,
%velocity relaxation time, destination position, type of
%destination (food=0, table=1, reached the table = 2)
attributes_p = 10; %[x y vx vy fx fy T_p obj_x obj_y d]
person = zeros(N_p,attributes_p);

% Relaxation time (for person-object interaction)
person(:,7) = 0.3*rand(size(person(:,7)))+0.5;

%Table attributes: position, table-to-people interaction constant, free or
%not (np=0 means free table, np%%, means that there are %% people at the
%table)
attributes_t = 3; %[x y np]
table = zeros(N_t,attributes_t);

%Table repulsion constant, minimum distance from person-table,maximum people at a table
table_const = [C_t min_dist_table max_p];

%Food attributes: position
attributes_f = 2; %[x y]
food = zeros(N_f,attributes_f);

% Definition of the (x,y) coordinate system (origin on bottom left)
% EMILE: the map can be change but remember that it is connected to the
% wall force field. If you change the limits or number of points, you also
% need to update the force field.
x_map = linspace(0,10,200);
y_map = linspace(0,10,200);
% Building the map
[X_map,Y_map] = meshgrid(x_map,y_map(end:-1:1));

% Initialization of the attributes of people, food, table and map
%People's initial position
min_distance = abs(X_map(1,1)-X_map(1,2))*10/N_p; % Minimum initial distance among people
[person(:,1), person(:,2)] = people_initial_position(person,X_map,Y_map,min_distance);
%People's initial velocity
%????
%People's initial destination is already set to zero (i.e the food)
%Peoples's relaxtion time
%Food location
food = [3.5 0.2; 6 0.2];
%Tables location
[table_x , table_y] = tablePositions(tableShape,N_t);
table(:,1:2) = [table_x table_y];

%The wall position is defined in the "wall_Emile script"

% Force field determined by the walls
[static_fx,static_fy] = readForceFile();
coloringTheMap(static_fx,static_fy,max(max(X_map))/1000,X_map,Y_map,1)
functionToPlotTheStaticField(X_map,Y_map,static_fx,static_fy)
%Defining the time steps
dt = 0.4;
final_time = 100;
% defining cost variables
s=0;
cost_v_tot=zeros(final_time/dt,N_p); %total cost velocity
cost_f_tot=zeros(final_time/dt,N_p); %total cost forces
cost_t1_tot=zeros(1,N_p); %total cost time to reach food
cost_t2_tot=zeros(1,N_p); %total cost time to reach tables

%MAIN LOOP:
for t=dt:dt:final_time
    % Initializing the forces to zero
    person(:,5:6) = zeros(size(person(:,5:6)));
    
    [x_obj, y_obj , table, person] = objective_direction(person,food,table,tableCapacity);
    person(:,8:9) = [x_obj y_obj]; % filling the destionations into the person matrix
    plotting_tables_food_people(food,table,person,X_map,Y_map,static_fx,static_fy,min_dist_table)
    
    %fromPositionToDirections(person,table,food);
    
    
    for i = 1:N_p
        
        people = person;
        people(i,:) = []; % getting rid of the row that belongs to i th person. this matrix is the matrix for the people other than i th person
        
        % Calculation of the forces. Previous verison looks OK. But it yields wrong answer. When it is used as A(x,y) = MultiOutFUnc(in)
        % only the first value of the output is used. [x,y]=MultiOutFunc
        % should be used
        
        [fx_wall ,fy_wall] = person_walls_force(person(i,:),static_fx,static_fy,X_map,Y_map,C_w);
        person(i,5:6) = person(i,5:6) + [fx_wall fy_wall];
        
        
        [fx_objective , fy_objective] = person_objective_force(person(i,:),v0,table,food);
        person(i,5:6) = person(i,5:6) + [fx_objective fy_objective];
        
        [fx_table, fy_table] = person_tables_force(person(i,:),table,food,table_const);
        person(i,5:6) = person(i,5:6) + [fx_table fy_table];
        
        [fx_people , fy_people] = person_people_force(person(i,:), people ,dt,A,B,sightAngle,sightCoef);
        person(i,5:6) = person(i,5:6) + [fx_people fy_people];
                
        [cost_t1,cost_t2] = cost_function_t1(i,person,cost_t1_tot,cost_t2_tot,t);
        cost_t1_tot(1,i)=cost_t1_tot(1,i)+cost_t1(1,i);
        cost_t2_tot(1,i)=cost_t2_tot(1,i)+cost_t2(1,i);
        if cost_t2_tot(1,i)~=0 && cost_t2_tot(1,i)>cost_t1_tot(1,i)
            cost_t2_tot(1,i)=cost_t2_tot(1,i)-cost_t1_tot(1,i);
        end
    end
    % Calculation of the position and velocity
    % Previous version give error because when it used as A(l,2:3) = func(), func gives only the first output
    [x_upt , y_upt , vx_upt , vy_upt] = update_position_velocity(person, dt,v_lim); % dummy variables
    person(:, 1:4) = [x_upt y_upt vx_upt vy_upt];
    [x_upt , y_upt , vx_upt , vy_upt, d_upt] = updating_objective(person,food,table,min_dist_obj,min_dist_table);
    person(:,[1:4 10]) = [x_upt , y_upt , vx_upt , vy_upt, d_upt];
    
    
    [cost_v,cost_f] = cost_function_v_f(person,vx_upt,vy_upt,v_lim,N_p);
    s=s+1;
    cost_v_tot(s,:)=cost_v; 
    cost_f_tot(s,:)=cost_f;
    
    pause(0.1)
end

cost_v_mean=sum(cost_v_tot)./sum(cost_v_tot~=0); %calculating cost of velocity for every person
cost_f_mean=sum(cost_f_tot)./sum(cost_v_tot~=0); %calculating cost of force for every person
cost_v_mean_max=max(cost_v_mean); %calculating max of velocity
cost_f_mean_max=max(cost_f_mean); %calculating max of force
cost_t1_tot_max=max(cost_t1_tot); %calculating max of cost to reach food table
cost_t2_tot_max=max(cost_t2_tot); %calculating max of cost to reach tables
%normalizing values in base 10
for i=1:N_p 
    cost_v_mean(i)=(cost_v_mean(i).*10)./(cost_v_mean_max);
    cost_f_mean(i)=(cost_f_mean(i).*10)./(cost_f_mean_max);
    cost_t1_tot(i)=(cost_t1_tot(i).*10)./(cost_t1_tot_max);
    cost_t2_tot(i)=(cost_t2_tot(i).*10)./(cost_t2_tot_max);
end
cost_total=round(sum(mean(cost_v_mean)+mean(cost_f_mean)+mean(cost_t1_tot)+mean(cost_t2_tot)),2); %calculating total cost

function [x , y] = tablePositions(tableShape,N_t)

x = zeros(N_t,1);
y = zeros(N_t,1);

if(tableShape == 1) % if in circular shape
    
    for i = 1:1:N_t
        
        x(i) = 4.75 + 1.5*cos(2*pi/N_t * (i-1));
        y(i) = 4 + 1.5*sin(2*pi/N_t * (i-1));
           
    end
    


elseif(tableShape ==2) % if it is rectangular
    
    x(1) = 4; x(2) = 4; x(3) = 6; x(4) = 6;
    y(1) = 6; y(2) = 2; y(3) = 2; y(4) = 6;
    
    for i = 1:1:(N_t-4)/2
       
        x(4+i) = 4;
        y(4+i) = 2 + 4/( (N_t-4)/2 + 1 )*i;
        
    end
    
    for i = 1:1:(N_t-4)/2
       
        x(N_t/2 + 2 + i) = 6;
        y(N_t/2 + 2 + i) = 2 + 4/( (N_t-4)/2 + 1 )*i;
        
    end
   
    
end    

end

function [x, y, vx, vy, d] = updating_objective(person,food,table,min_dist_obj,min_dist_table)
N_t = size(table,1);
N_f = size(food,1);
N_p = size(person,1);
%Reasoning: the position of the people has been updated. If the person if
%too close to the objective, it is stopped, positioned at a safety distance
%and the objective is changed (from d=0 food to d=1 table or d=???? table
%reached)
corner = [2.5 7; 7 7];
dist_person_obj = ((person(:,1)-person(:,8)).^2+(person(:,2)-person(:,9)).^2).^(1/2);
dist_person_obj( (person(:,8)==corner(1,1) & person(:,9)==corner(1,2)) | (person(:,8)==corner(2,1) & person(:,9)==corner(2,2)) ) = max(dist_person_obj);
min_dist_obj = repelem(min_dist_obj,N_p,1) + (person(:,10)~=0).*min_dist_table;

v_tot = (person(:,3).^2+person(:,4).^2).^(1/2);
%If the person has reached the objective and is too near to it, he/she is
%shift to the min_dist from the objective

x = (dist_person_obj>min_dist_obj).*person(:,1)+(dist_person_obj<=min_dist_obj).*(person(:,8)-min_dist_obj.*person(:,3)./v_tot);
y = (dist_person_obj>min_dist_obj).*person(:,2)+(dist_person_obj<=min_dist_obj).*(person(:,9)-min_dist_obj.*person(:,4)./v_tot);

%If the person has reached the objective, he/she stops.
vx = (dist_person_obj>min_dist_obj).*person(:,3);
vy = (dist_person_obj>min_dist_obj).*person(:,4);

%If the person is too close to the food (d=0), the objective becomes the
%table (d = d+1 = 1)

d = person(:,10) + (dist_person_obj<=min_dist_obj & (person(:,10) == 0))*100; % becomes 100 when the food is first taken and then it will be the table number

%d = (dist_person_obj<=min_dist_obj & (person(:,10) == 0))*sth % write
%what you want d to be after people reach

end

function [fx_table, fy_table] = person_tables_force(person,table,food,table_const)
%This functions calculates the repulsion force of the table 

% Number of tables
N_t = size(table,1);

% Excluding the table that is the objective of one person (if he/she has
% already taken the food.

[dx, dy, nearest_objective] = fromPositionToDirections(person,table,food);
if person(10)~=0
    table = table(1:N_t~=nearest_objective,:);
end

dist_people_table = sqrt( (person(1)-table(:,1)).^2 + (person(2)-table(:,2)).^2);
fx_table = sum(table_const(1).*(person(1)-table(:,1))./(dist_people_table.^(4/2)));
fy_table = sum(table_const(1).*(person(2)-table(:,2))./(dist_people_table.^(4/2)));

end

function [x_obj, y_obj, table, person] = objective_direction(person,food,table,tableCapacity)
%For now, the direction of the objective is simply the direction towards the nearest food
% (or nearest table, according to people's d, the objective).
%Neglecting the full table
%table = table(table(:,3)==0,:);

N_t = size(table(:,1),1);
N_f = size(food(:,1),1);
%Number of people that want to go toward a table of food
N_p_food = size(person(person(:,10)==0),1);
N_p_table = size(person(person(:,10)~=0),1);
% Division of the people that want to go toward a table of food into two
%matrices
person_food = person(person(:,10)==0,:);
person_table = person(person(:,10)~=0,:);

%
%Every row is a person, every column a food point / table
dist_food_person = sqrt((repelem(person_food(:,1),1,N_f)-repelem(food(:,1)',N_p_food,1)).^2+(repelem(person_food(:,2),1,N_f)-repelem(food(:,2)',N_p_food,1)).^2);
dist_table_person = sqrt((repelem(person_table(:,1),1,N_t)-repelem(table(:,1)',N_p_table,1)).^2+(repelem(person_table(:,2),1,N_t)-repelem(table(:,2)',N_p_table,1)).^2);

%If you are exactly on a table of on a food point, a numerical error would
%appear. In order to avoid this problem, zero distances are approximanted
%to a very small number
dist_food_person(dist_food_person==0) = dist_food_person(dist_food_person==0)+10*1e-10;
dist_table_person(dist_table_person==0) = dist_table_person(dist_table_person==0)+10*1e-10;


%Identification of the nearest food point to each person
nearest_food = (dist_food_person'== min(dist_food_person'))'*(1:N_f)';
%Identification of the objective table to each person
%nearest_table = (dist_table_person'==min(dist_table_person'))'*(1:N_t)';

[sortedTableDist , sortedTableIndex] = sort(dist_table_person,2);

table_person_indices = find(person(:,10)); % indice number of poeple whose objective is table

for i = 1:1:size(table_person_indices,1) % iterating through table_people
    for j = 1:1:N_t % iterating through tables
        
    if( table( sortedTableIndex(i,j) , 3)  < tableCapacity && person(table_person_indices(i),10) == 100 ) % if closest table is not full
        
        table( sortedTableIndex(i,j) , 3) = table( sortedTableIndex(i,j) , 3) + 1; % increasing the value indicating how many people on that table
        
        person(table_person_indices(i),10) = sortedTableIndex(i,j);
        
        break;
    end
    
    end
    
end

x_obj(person(:,10)==0) = food(nearest_food,1);
y_obj(person(:,10)==0) = food(nearest_food,2);

x_obj(person(:,10)~=0 ) = table(person(table_person_indices,10),1);
y_obj(person(:,10)~=0 ) = table(person(table_person_indices,10),2);

% what happens
% x_obj(person(:,10)==condition that table is reached) = objective after table is reahed ??;
% y_obj(person(:,10)==condition that table is reached) = objective after table is reahed ??;


%We prefer colunm vertors
x_obj = reshape(x_obj,[],1); 
y_obj = reshape(y_obj,[],1);



%If a person is situated in one of the top corner, the the destination is
%splitted into a polygonal line going to the corner, then to food
corner = [2.5 7; 7 7];

%If you are in the bad corner of the room but you can see the second food,
%move towards it.
x_obj((person(:,2)>=corner(1,2) & person(:,2)<food(1,2)+(person(:,1)-food(1,1)).*(corner(1,2)-food(1,2))/(corner(1,1)-food(1,1)) & person(:,2)>=food(2,2)+(person(:,1)-food(2,1)).*(corner(1,2)-food(2,2))/(corner(1,1)-food(2,1)))) = food(2,1);
y_obj((person(:,2)>=corner(1,2) & person(:,2)<food(1,2)+(person(:,1)-food(1,1)).*(corner(1,2)-food(1,2))/(corner(1,1)-food(1,1)) & person(:,2)>=food(2,2)+(person(:,1)-food(2,1)).*(corner(1,2)-food(2,2))/(corner(1,1)-food(2,1)))) = food(2,2);
x_obj((person(:,2)>=corner(2,2) & person(:,2)<food(2,2)+(person(:,1)-food(2,1)).*(corner(2,2)-food(2,2))/(corner(2,1)-food(2,1)) & person(:,2)>=food(1,2)+(person(:,1)-food(1,1)).*(corner(2,2)-food(1,2))/(corner(2,1)-food(1,1)))) = food(1,1);
y_obj((person(:,2)>=corner(2,2) & person(:,2)<food(2,2)+(person(:,1)-food(2,1)).*(corner(2,2)-food(2,2))/(corner(2,1)-food(2,1)) & person(:,2)>=food(1,2)+(person(:,1)-food(1,1)).*(corner(2,2)-food(1,2))/(corner(2,1)-food(1,1)))) = food(1,2);

%If you cannot see nono of the food points, go towards the corner of the
%room in order to minimize the the path
x_obj((person(:,2)>=corner(1,2) & person(:,2)<food(2,2)+(person(:,1)-food(2,1)).*(corner(1,2)-food(2,2))/(corner(1,1)-food(2,1)))) = corner(1,1);
y_obj((person(:,2)>=corner(1,2) & person(:,2)<food(2,2)+(person(:,1)-food(2,1)).*(corner(1,2)-food(2,2))/(corner(1,1)-food(2,1)))) = corner(1,2);
x_obj((person(:,2)>=corner(2,2) & person(:,2)<food(1,2)+(person(:,1)-food(1,1)).*(corner(2,2)-food(1,2))/(corner(2,1)-food(1,1)))) = corner(2,1);
y_obj((person(:,2)>=corner(2,2) & person(:,2)<food(1,2)+(person(:,1)-food(1,1)).*(corner(2,2)-food(1,2))/(corner(2,1)-food(1,1)))) = corner(2,2);




end

function [fx, fy] = person_objective_force(person,v0,table,food)
%This function is to calculate the forces between
%person and tables, person and foods. We assume that people are not only
%attracted by nearest foods/tables, but also attracted by other relatively
%distant foods and tables, with which it is possible for people to change
%the directions of their heading because of the interference of other
%people.

% %Neglecting the full table
% table = table(table(:,4)==0,:);
% %Number of tables and foods we should consider
% N_t = max(size(table(:,1)));
% N_f = max(size(food(:,1)));
% 
% % Division of the people that want to go toward a table of food into two
% %matrices
% person_food = person(person(10)==0,:);
% person_table = person(person(10)==1,:);
% 
% %Calculate the distance between person and tables/foods
% dist_person_food_x = repelem(person_food(1),1,N_f)-food(:,1)';
% dist_person_food_y = repelem(person_food(2),1,N_f)-food(:,2)';
% dist_person_food = sqrt( dist_person_food_x.^2 + dist_person_food_y.^2 );

% dist_person_table_x = repelem(person_table(1),1,N_t)-table(:,1)'; % Nf are changed with Nt
% dist_person_table_y = repelem(person_table(2),1,N_t)-table(:,2)'; % This is correct right ??

% if person(10)==0 %return force caused by foods
%     C_f = 0.01; %Food to person force rate
%     fx = -C_f*abs(sum(dist_person_food_x./(abs(dist_person_food)).^3)); 
%     fy = -C_f*abs(sum(dist_person_food_y./(abs(dist_person_food)).^3));
%     
%     
% elseif person(10)==1 %return force caused by tables . % SHOULD BE ALSO CORRECTED. THIS IS WRONG Hangxi: change it into -C_t
%     C_t = 0.01*max(table(:, 3));
%     fx = -C_t*abs(sum(dist_person_table_x./(abs(dist_person_table_x)).^3));
%     fy = -C_t*abs(sum(dist_person_table_y./(abs(dist_person_table_y)).^3));
%end
[dx, dy] = fromPositionToDirections(person,table,food);

T_relaxation = person(7);
fx = (v0*dx-person(3))/T_relaxation;
fy = (v0*dy-person(4))/T_relaxation;
end

function [fx, fy] = person_people_force(person, people,dt,A,B,sightAngle,sightCoef)
%This function is to calculate the force caused by all other people to person A
%Following is the distance between two people, every element represents the
%distance between the person and the other people


N_p = size(people(:,1),1);
dist_person_people_x = (repelem(person(1),1,N_p)-people(:,1)')';
dist_person_people_y = (repelem(person(2),1,N_p)-people(:,2)')';

dist_person_people = sqrt(dist_person_people_x.^2 + dist_person_people_y.^2); % total distance between 2 persons

y_x = (people(:,3) - person(3)*ones(size(people,1),1))*dt;
y_y = (people(:,4) - person(4)*ones(size(people,1),1))*dt;

b = sqrt( ( dist_person_people + sqrt( (dist_person_people_x - (people(:,3) - person(3))*dt).^2 + (dist_person_people_y - (people(:,4) - person(4))*dt).^2 ) ).^2 - (y_x.^2 + y_y.^2) ) / 2;

f_part1 = A * exp(-b/B); % The repulsive potential
f_part2 = (dist_person_people + sqrt((dist_person_people_x - y_x).^2 + (dist_person_people_y - y_y).^2 ) ) ./ (2*b);
f_part3x = (dist_person_people_x ./ dist_person_people) + ( (dist_person_people_x - y_x) ./ sqrt( (dist_person_people_x - y_x).^2 + (dist_person_people_y - y_y).^2 ) );
f_part3y = (dist_person_people_y ./ dist_person_people) + ( (dist_person_people_y - y_y) ./ sqrt( (dist_person_people_x - y_x).^2 + (dist_person_people_y - y_y).^2 ) );

% arranging the force componenets. such that it has the direction of the distance between people


coef = double(abs(person_people_angle(person,people)) < sightAngle);
   
coef(coef == 0) = sightCoef; % if it is out of sight, coef becomes sight coef
 
fx = coef.*(f_part1 .* f_part2 * 0.5 .* f_part3x); % array that shows the force from other people
fy = coef.*(f_part1 .* f_part2 * 0.5 .* f_part3y);

fx = sum(fx); % sum of those forces
fy = sum(fy);

end

function gama = person_people_angle(person,people)

% This function is to calculate the angles between person and all other people

velocityAngle = atan2(person(4) , person(3));

gamaInertial = atan2(people(:,2) - person(2)*ones(size(people,1),1) , people(:,1) - person(1)*ones(size(people,1),1)  ); % angle between person and other poeple in inertail fram

gama = gamaInertial - velocityAngle;


end

function [x, y, vx, vy] = update_position_velocity(person, dt , v_lim)

vx = person(:,3) + dt*person(:,5); % x velocity of a person
vy = person(:,4) + dt*person(:,6); % y veloctiy of a person
v = sqrt(vx.^2 + vy.^2); % speed  of a person 



vx(v>v_lim) = v_lim .* vx(v>v_lim) ./ v(v>v_lim);
vy(v>v_lim) = v_lim .* vy(v>v_lim) ./ v(v>v_lim);

x = person(:,1) + dt*vx;
y = person(:,2) + dt*vy;

end

function [fx,fy]=person_walls_force(person,static_fx,static_fy,X_map,Y_map,C_w)

% The force applied by wall is an interpolation among the forces of the
% nearest four points of the static field. In order to find four points, I
% want that the position of each person doens't coincide with the edges of
% the sides of the square connected the four points.
%Approximation: if this condition is not satisfied, I slightly move the
%person (in the following algorithm, I move it on top left).
bolean_variable = 1;

while bolean_variable == 1
    if sum(person(1)==X_map)~=0
        person(1) = person(1) + max(max(X_map))/1e6;
    elseif sum(person(2)==Y_map)~=0
        person(2) = person(2) + max(max(Y_map))/1e6;
    else
        bolean_variable = 0;
    end
end

% Now I need to find the forces on the four edges.

%Distance between points in the map
dist_x = abs(X_map(1,2)-X_map(1,1));
dist_y = abs(Y_map(2,1)-Y_map(1,1));

if dist_x==0 || dist_y ==0
    disp('problem in calculating the distances in the function person_wallInteraction')
end

% Looking for the position of the four exges
pos_x = (X_map<person(1)+dist_x & X_map>person(1)-dist_x);
pos_y = (Y_map<person(2)+dist_y & Y_map>person(2)-dist_y);
pos_xy =  logical(pos_x.*pos_y);
% Position of the four edges
x_edges = X_map(pos_xy);
y_edges = Y_map(pos_xy);
% Forces at the two edges
fx_edges = static_fx(pos_xy);
fy_edges = static_fy(pos_xy);


distPerson_edge = ( (x_edges-person(1)).^2   + (y_edges-person(2)).^2 ).^(1/2);
% Calculating the average force applied on the person located among the
% four edges. The force is calculated as a weigthed average based on the
% distance of the person from the four points.
fx = C_w*sum(fx_edges.*distPerson_edge)/sum(distPerson_edge);
fy = C_w*sum(fy_edges.*distPerson_edge)/sum(distPerson_edge);
end

function [static_fx,static_fy] = readForceFile()
fileID = fopen('static_fx.txt','r');
static_fx =fscanf(fileID,'%f');
static_fx = reshape(static_fx,[sqrt(max(size(static_fx))),sqrt(max(size(static_fx)))]);
fclose(fileID);

fileID = fopen('static_fy.txt','r');
static_fy =fscanf(fileID,'%f');
static_fy = reshape(static_fy,[sqrt(max(size(static_fy))),sqrt(max(size(static_fy)))]);
fclose(fileID);
end

function coloringTheMap(static_fx,static_fy,distanceAmongPointsOnTheWall,X_map,Y_map,wallRepulsionConstant)

statif_tot_f = (static_fx.^2+static_fy.^2).^(1/2);
Fmax = wallRepulsionConstant/(distanceAmongPointsOnTheWall^2);
Fmin = wallRepulsionConstant/(max(max(X_map)/2)^2);
%Smoothing the force field
statif_tot_f(statif_tot_f>Fmax) = Fmax;
statif_tot_f(statif_tot_f<Fmin) = Fmin;

figure(1)
hold on
%Coloring the map
% Create a sample image of a ramp.
% Smoothing the map
heatMap = (statif_tot_f./Fmax).^(1/5);


% Display it.
imagesc(X_map(1,1:end),Y_map(1:end,1),heatMap);
%{
% Initialize a color map array of 101 colors.
colorMap = jet(30);
% Apply the colormap and show the colorbar
colormap(colorMap);
%}
colorbar;
axis image;
hold off
end

function functionToPlotTheStaticField(X_map,Y_map,static_fx,static_fy)
figure(1)
hold on
%Transformation of forces on unit verctors
fx = static_fx./sqrt(static_fx.^2+static_fy.^2);
fy = static_fy./sqrt(static_fx.^2+static_fy.^2);
quiver(X_map(1:3:end,1:3:end),Y_map(1:3:end,1:3:end),fx(1:3:end,1:3:end),fy(1:3:end,1:3:end),0.6)
hold off
end

function plotting_tables_food_people(food,table,person,X_map,Y_map,static_fx,static_fy,min_dist_table)

[dx, dy] = fromPositionToDirections(person,table,food);

angle = linspace(0,2*pi,20);
x_around_table = min_dist_table.*cos(angle);
y_around_table = min_dist_table.*sin(angle);

figure(1)
coloringTheMap(static_fx,static_fy,max(max(X_map))/1000,X_map,Y_map,1)
title('Plot of the direction of the moving force (at first it points toward the food)')
hold on

plot(food(:,1),food(:,2),'r+',table(:,1),table(:,2),'ro',person(:,1),person(:,2),'g*')
%{
plot(food(:,1),food(:,2),'r+',person(:,1),person(:,2),'g*')
for i = 1:size(table,1)
    plot((x_around_table+table(i,1)),(y_around_table+table(i,2)),'r');
end
%}

%dx(person(:,10)==condition for people reach table) = 0;
%dy(person(:,10)==condition for people reach table) = 0;
quiver(person(:,1),person(:,2),dx,dy,0.4,'r')
axis([min(min(X_map)) max(max(X_map)) min(min(Y_map)) max(max(Y_map))])

hold off
end

function [pos_x, pos_y] = people_initial_position(person,X_map,Y_map,min_distance)
N_p = size(person,1);
boolean_variable = 1;
while boolean_variable==1
    person(1:floor(N_p/2),1) = 2*rand(size(person(1:floor(N_p/2),1)))+0.5;
    person(floor(N_p/2)+1:end,1) = 2*rand(size(person(floor(N_p/2)+1:end,1)))+7.5;
    person(:,2) = 2*rand(size(person(:,2)))+7.5;
    dist_person_person_x_map = (repelem(person(:,1),1,N_p)-repelem(person(:,1),1,N_p)')';
    dist_person_person_y_map = (repelem(person(:,2),1,N_p)-repelem(person(:,2),1,N_p)')';
    dist_person_person_x = dist_person_person_x_map(dist_person_person_x_map~=0 & dist_person_person_y_map~=0)';
    dist_person_person_x = reshape(dist_person_person_x,[],N_p)';
    dist_person_person_y = dist_person_person_y_map(dist_person_person_x_map~=0 & dist_person_person_y_map~=0)';
    dist_person_person_y = reshape(dist_person_person_y,[],N_p)';
    tot_dist = (dist_person_person_x.^2+dist_person_person_y.^2).^(3/2);
    if sum(tot_dist<min_distance)==0
        boolean_variable=0;
    end
end
pos_x = person(:,1);
pos_y = person(:,2);
end

function [dx, dy, nearest_objective] = fromPositionToDirections(person,table,food)
dist_person_objective_x = person(:,8) - person(:,1);
dist_person_objective_y = person(:,9) - person(:,2);
tot_dist_person_objective = (dist_person_objective_x.^2+dist_person_objective_y.^2).^(1/2);
dx = dist_person_objective_x./tot_dist_person_objective;
dy = dist_person_objective_y./tot_dist_person_objective;

%Dividing the people into the ones looking for food or tables
person_food = person(person(:,10)==0,:);
person_table = person(person(:,10)~=0,:);

N_t = size(table,1);
N_f = size(food,1);
N_p_food = size(person_food,1);
N_p_table = size(person_table,1);
%Finding the index of the objective (to recognize it this index represent
%food or table, look at person(:,10)  )

nearest_objective_food = (repelem(person_food(:,8),1,N_f) == repelem(food(:,1)',N_p_food,1) & repelem(person_food(:,9),1,N_f) == repelem(food(:,2)',N_p_food,1))*(1:N_f)';
nearest_objective_table = (repelem(person_table(:,8),1,N_t) == repelem(table(:,1)',N_p_table,1) & repelem(person_table(:,9),1,N_t) == repelem(table(:,2)',N_p_table,1))*(1:N_t)';

nearest_objective(person(:,10)==0) = nearest_objective_food;
nearest_objective(person(:,10)~=0) = nearest_objective_table;

end

%the next two functions are cost functions. The cost of each different
%combination of parameters can be determined by considering 4 things: the
%velocity of each person,the force applied to each person, the time required to reach the food table and the
%time required to find an empty table.
function [cost_v,cost_f] = cost_function_v_f(person,vx_upt,vy_upt,v_lim,N_p)
%function to calculate the cost considering the velocity of each person and the force 
% applied to each person. 
% The cost for a general velocity is set to be proportional to 1/v ( as used often in
% literature). When the velocity of a person is equal to v_lim then the
% cost is equal to 0. The function returns an array of values with the cost
% of every person considering the velocity for every dt and an array of
% values of every person considering the force applied every dt.
v=sqrt(vx_upt.^2+vy_upt.^2); %velocity
cost_v=1./v;
for p=1:1:N_p
    if v(p)==v_lim
        cost_v(p)=0;
    end
    if v(p)==0
        cost_v(p)=1/0.1;
    end
end
f=sqrt((person(:,5)).^2+(person(:,6)).^2); %force
cost_f=f;
end

function [cost_t1,cost_t2] = cost_function_t1(i,person,cost_t1_tot,cost_t2_tot,t)
% function to calculate the cost considering the time to reach the
% objectives of every person.One fist cost consideres the time to reach the
% food table and the cost is equal to the time to reach it. The second cost
% considers the time to find and reach an empty table. 
if person(i,10)==1 && cost_t1_tot(1,i)==0
    time1(i)=t; 
    cost_t1(1,i)=t;
else 
    cost_t1(1,i)=0;
end

if person(i,10)==2 && cost_t2_tot(1,i)==0
    time2(i)=t;
    cost_t2(1,i)=t;
else
    cost_t2(1,i)=0;
   
end 
end
