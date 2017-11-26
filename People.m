close all
clear
clc
% this file is a sample program to illustract the structure of our social
% model algorithm for Apero

N_p = 5;   % number of people
N_t = 10;  % number of table
N_f = 2;   % number of foods

%People's attributes: position, velocity, undergoing force,
%velocity relaxation time, destination direction, type of
%destination (food=0, table=1).
attributes_p = 10; %[x y vx vy fx fy T_p dx dy d]
person = zeros(N_p,attributes_p);

%Table attributes: position, table-to-people interaction constant, free or
%not (free = 0, full =1)
attributes_t = 4; %[x y C_t f]
table = zeros(N_t,attributes_t);

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
person(1:floor(N_p/2),1) = 2.5*rand(size(person(1:floor(N_p/2),1)));
person(floor(N_p/2)+1:end,1) = 2.5*rand(size(person(floor(N_p/2)+1:end,1)))+7.5;
person(:,2) = 3*rand(size(person(:,2)))+7;
%People's initial velocity
%????
%People's initial destination is already set to zero (i.e the food)
%Peoples's relaxtion time
person(:,7) = 5*rand(size(person(:,7)));
%Food location
food = [4 0; 6 0];
%Tables location
table(:,1:2) = 3*rand(size(table(:,1:2)))+3;
%At first, all the tables are free
%The wall position is defined in the "wall_Emile script"

% Force field determined by the walls
[static_fx,static_fy] = readForceFile();
coloringTheMap(static_fx,static_fy,max(max(X_map))/1000,X_map,Y_map,1)
functionToPlotTheStaticField(X_map,Y_map,static_fx,static_fy)
%Defining the time steps
dt = 0.4;
final_time = 50;

%MAIN LOOP:
for t=dt:dt:final_time
    % Initializing the forces to zero
    person(:,5:6) = zeros(size(person(:,5:6)));
    
    % EMILE: "In my opinion the people loop can be avoided, anyway to start
    % with a simple code probably it's better to use it. If you are able,
    % try to build the functions so that they can calculate the outcomes
    % also if they receive as input all the people at the same time. I did
    % it for the objective research, it was fun but hard.Good luck!"
    [dx, dy] = objective_direction(person,food,table);
    person(:,8:9) = [dx' dy']; % filling the destionations into the person matrix
    
    for i = 1:N_p
        
        people = person;
        people(i,:) = []; % getting rid of the row that belongs to i th person. this matrix is the matrix for the people other than i th person
        
        % Calculation of the forces. Previous verison looks OK. But it yields wrong answer. When it is used as A(x,y) = MultiOutFUnc(in)
        % only the first value of the output is used. [x,y]=MultiOutFunc
        % should be used
        
        %[fx_wall ,fy_wall] = person_walls_force(person(i,:),static_fx,static_fy,X_map,Y_map)
        %person(i,5:6) = person(i,5:6) + [fx_wall fy_wall];
        
        [fx_objective , fy_objective] = person_objective_force(person(i,:));
        person(i,5:6) = person(i,5:6) + [fx_objective fy_objective];
        
        %[fx_table fy_table] = person_tables_force(person(i,:),other_stuff);
        %person(i,5:6) = person(i,5:6) + [fx_table fy_table];
        
        [fx_people , fy_people] = person_people_force(person(i,:), people ,dt);
        person(i,5:6) = person(i,5:6) + [fx_people fy_people];
        
    end
    % Calculation of the position and velocity
    % Previous version give error because when it used as A(l,2:3) = func(), func gives only the first output
    [x_upt , y_upt , vx_upt , vy_upt] = update_position_velocity(person, dt); % dummy variables
    person(:, 1:4) = [x_upt y_upt vx_upt vy_upt];
    
    pause(1)
    
end



function [dx, dy] = objective_direction(person,food,table)
%For now, the direction of the objective is simply the direction towards the nearest food
% (or nearest table, according to people's d, the objective).
%Neglecting the full table
table = table(table(:,4)==0,:);

N_t = max(size(table(:,1)));
N_f = max(size(food(:,1)));
%Number of people that want to go toward a table of food
N_p_food = max(size(person(person(:,10)==0)));
N_p_table = max(size(person(person(:,10)==1)));
% Division of the people that want to go toward a table of food into two
%matrices
person_food = person(person(:,10)==0,:);
person_table = person(person(:,10)==1,:);

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
nearest_food = (dist_food_person'==min(dist_food_person'))'*(1:N_f)';
%Identification of the nearest table to each person
nearest_table = (dist_table_person'==min(dist_table_person'))'*(1:N_t)';

% Unit vectors pointing towards the nearest food
dx_food = (food(nearest_food,1)-person_food(:,1))./min(dist_food_person')';
dy_food = (food(nearest_food,2)-person_food(:,2))./min(dist_food_person')';
% Unit vectors pointing towards the nearest table
dx_table = (table(nearest_table,1)-person_table(:,1))./min(dist_table_person')';
dy_table = (table(nearest_table,2)-person_table(:,2))./min(dist_table_person')';

%disp('person(:,10)==0 means food, person(:,10)==1 means table')
dx(person(:,10)==0) = dx_food;
dx(person(:,10)==1) = dx_table;
dy(person(:,10)==0) = dy_food;
dy(person(:,10)==1) = dy_table;

figure(1)
title('Plot of the direction of the moving force (at first it points toward the food)')
hold on
plot(food(:,1),food(:,2),'r+',table(:,1),table(:,2),'ro',person(:,1),person(:,2),'g*')
quiver(person(:,1),person(:,2),dx',dy',0.4)
hold off
end

function [fx, fy] = person_objective_force(person)
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

T_relaxation = person(7);
v0 = 1; 
fx = (v0*person(8)-person(3))/T_relaxation;
fy = (v0*person(9)-person(4))/T_relaxation;
end

function [fx, fy] = person_people_force(person, people,dt)
%This function is to calculate the force caused by person B to person A
%Following is the distance between two people, every element represents the
%distance between the person and the other people

A = 1; % Taken same as the paper
B = 0.2; % same as paper

N_p = max(size(people(:,1)));
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
fx = sum(f_part1 .* f_part2 * 0.5 .* f_part3x);
fy = sum(f_part1 .* f_part2 * 0.5 .* f_part3y);

end

function [x, y, vx, vy] = update_position_velocity(person, dt)
vx = person(:,3) + dt*person(:,5);
vy = person(:,4) + dt*person(:,6);
x = person(:,1) + dt*vx;
y = person(:,2) + dt*vy;
end

function [fx,fy]=person_walls_force(person,static_fx,static_fy,X_map,Y_map)

% The force applied by wall is an interpolation among the forces of the
% nearest four points of the static field. In order to find four points, I
% want that the position of each person doens't coincide with the edges of
% the sides of the square connected the four points.
%Approximation: if this condition is not satisfied, I slightly move the
%person (in the following algorithm, I move it on top left).
bolean_variable = 1;

coef = 0.07; % this coefficient is added, because wall forces are huuge.
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
fx = coef*sum(fx_edges.*distPerson_edge)/sum(distPerson_edge);
fy = coef*sum(fy_edges.*distPerson_edge)/sum(distPerson_edge);
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
