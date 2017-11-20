% this file is a sample program to illustract the structure of our social
% model algorithm for Apero

N_p = 20;   % number of people
N_t = 10;  % number of table
N_f = 2;   % number of foods

%People's attributes: position, velocity, undergoing force,
%destination direction, person-to-person interaction constant.
attributes_p = 9; %[x y vx vy fx fy dx dy C]
person = zeros(N_p,attributes_p);

%Table attributes: position, table-to-people interaction constant
attributes_t = 3; %[x y C]
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
% Building the map (origin on bottom left)
[X_map,Y_map] = meshgrid(x_map,y_map(end:-1:1));

% Initialization of the attributes of people, food, table and map
%People's initial position
%????
%People's initial velocity
%????
%Food location
%????
%Table location
%????
%The wall position is defined in the "wall_Emile script"

% Force field determined by the walls
[static_fx,static_fy] = readForceFile();

%Defining the time steps
dt = 0.1;
final_time = 10;

%MAIN LOOP:
for t=dt:dt:final_time
    % Initializing the forces to zero
    person(:,5:6) = zeros(size(person(:,5:6)));
    
    % EMILE: "In my opinion this for loop can be avoided, anyway to start
    % with a simple code probably it's better to use it. If you are able,
    % try to build the functions so that they can calculate the outcomes
    % also if they receive as input all the people at the same time. Good
    % luck!"
    for i = 1:N_p
        % Updating the direction of the objective driving force
        % person(i,7:8) = objective_direction(person(i,:),other_stuff);
        
        % Calculation of the forces
        person(i,5:6) = person(i,5:6) + person_walls_force(person(i,:),static_fx,static_fy,X_map,Y_map);
        %person(i,5:6) = person(i,5:6) + person_objective_force(person(i,:),other_stuff);
        %person(i,5:6) = person(i,5:6) + person_tables_force(person(i,:),other_stuff);
         
        
        % Calculation of the position and velocity
        %person(i,1:4) = position_velocity(person(i,:),other_stuff);
    end
end



function [dx dy] = objective_direction(person,other_stuff)

end

function [fx fy] = person_objective_force(person,other_stuff)

end

function [fx fy] = person_tables_force(person,other_stuff)

end

function [x y vx vy] = position_velocity(person,other_stuff)

end

function [fx,fy]=person_walls_force(person,static_fx,static_fy,X_map,Y_map)

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
fx = sum(fx_edges.*distPerson_edge)/sum(distPerson_edge);
fy = sum(fy_edges.*distPerson_edge)/sum(distPerson_edge);
end
function [static_fx,static_fy] = readForceFile()
fileID = fopen('wall_fx.txt','r');
static_fx =fscanf(fileID,'%f');
static_fx = reshape(static_fx,[sqrt(max(size(static_fx))),sqrt(max(size(static_fx)))]);
fclose(fileID);

fileID = fopen('wall_fy.txt','r');
static_fy =fscanf(fileID,'%f');
static_fy = reshape(static_fy,[sqrt(max(size(static_fy))),sqrt(max(size(static_fy)))]);
fclose(fileID);
end
