clc
close all
clear

% Definition of the (x,y) coordinate system (on bottom left)
x_map = linspace(0,10,100);
y_map = linspace(10,0,100);

% Building the map
[X_map,Y_map] = meshgrid(x_map,y_map);

% Initializing the static force field to zero
static_fx = zeros(max(size(x_map)));
static_fy = zeros(max(size(x_map)));

% Position of the walls. The 2x2 matrix "wall" is made of two rows, each of
% those contains the (x,y) coordinates of the edges of the wall.
wall(1:2,1:2,1) = [2.5 0; 7 0];
wall(1:2,1:2,2) = [7 0; 7 7];
wall(1:2,1:2,3) = [7 7; 10 7];
wall(1:2,1:2,4) = [10 7; 10 10];
wall(1:2,1:2,5) = [10 10; 0 10];
wall(1:2,1:2,6) = [0 10; 0 7];
wall(1:2,1:2,7) = [0 7; 2.5 7];
wall(1:2,1:2,8) = [2.5 7; 2.5 0];



% Distance between points along the walls (constant for all walls)
distanceAmongPointsOnTheWall = max(max(X_map))/1000;
% Constant in front of the force 
disp('definition of the force can be changed')
wallRepulsionConstant = 1;

% Calculation of the static force gradient in the map determined by the
% repulsive forces of the walls
disp('It is better to save the static forces matrices in a file and just open the file in the main script')
for i=1:max(size(wall(1,1,:)))
    [static_fx,static_fy] = staticForceFunction(wall(:,:,i),X_map,Y_map,static_fx,static_fy,distanceAmongPointsOnTheWall,wallRepulsionConstant);
end
functionToPlotTheStaticField(X_map,Y_map,static_fx,static_fy)

% Example of a person in the map
person = [5 8];
% Calculation (and plot) of the forces applied by wall on that person
[wall_fx,wall_fy]=person_wallInteraction(person,static_fx,static_fy,X_map,Y_map);














%FUNCTIONS
function [static_fx,static_fy] = staticForceFunction(wall,X_map,Y_map,static_fx,static_fy,distanceAmongPointsOnTheWall,wallRepulsionConstant)

% Discretization of the wall
lengthOfTheWall = sqrt( (wall(1,1)-wall(2,1))^2 + (wall(1,2)-wall(2,2))^2);
numberOfPointsOnTheWall = floor(lengthOfTheWall/distanceAmongPointsOnTheWall);
angular_coeff = (wall(1,2)-wall(2,2))/(wall(1,1)-wall(2,1));
wall_x_end = wall(1,1)+ lengthOfTheWall/sqrt(1+angular_coeff^2)*(-1)^( wall(2,1)<wall(1,1));
wall_y_end = wall(1,2)+ lengthOfTheWall/sqrt(1+1/angular_coeff^2)*(-1)^( wall(2,2)<wall(1,2));

%Now the wall is discretized into several points. The distance between the
% points is constant for all walls.
x_wall = linspace(wall(1,1),wall_x_end,numberOfPointsOnTheWall);
y_wall = linspace(wall(1,2),wall_y_end,numberOfPointsOnTheWall);

%Plot of the discretized wall
figure(1)
hold on
plot(x_wall,y_wall,'r')
hold off

% Small approximation: no points of the wall can coincide with a point in 
%the map, otherwise I'll get a zero distance between point and wall, i.e.
%an infinite force.
% If it occurs, the the wall is slightly moved.
bolean_variable = 1;
while(bolean_variable ==1)
    for i = 1:numberOfPointsOnTheWall
        bolean_variable = 1;
        if sum(Y_map(X_map==x_wall(i))==wall(1,2))~=0
            x_wall = linspace(wall(1,1)+max(max(X_map))/1e6,wall_x_end+max(max(X_map))/1e6,numberOfPointsOnTheWall);
            continue;
        else
            bolean_variable = 0;
        end
    end
end


% Calculation of the effect of the wall
% Probably this part can be optimized by using the function gradient
for k=1:numberOfPointsOnTheWall
    static_fx = static_fx + wallRepulsionConstant.*(X_map-x_wall(k))./((X_map-x_wall(k)).^2+(Y_map-y_wall(k)).^2).^(3/2);
    static_fy  = static_fy + wallRepulsionConstant.*(Y_map-y_wall(k))./((X_map-x_wall(k)).^2+(Y_map-y_wall(k)).^2).^(3/2);
end



end

%Function that plots the wall
function functionToPlotTheStaticField(X_map,Y_map,static_fx,static_fy)
figure(1)
hold on
%Transformation of forces on unit verctors
fx = static_fx./sqrt(static_fx.^2+static_fy.^2);
fy = static_fy./sqrt(static_fx.^2+static_fy.^2);
quiver(X_map(1:3:end,1:3:end),Y_map(1:3:end,1:3:end),fx(1:3:end,1:3:end),fy(1:3:end,1:3:end),0.6)
hold off
end


function [wall_fx,wall_fy]=person_wallInteraction(person,static_fx,static_fy,X_map,Y_map)

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
wall_fx = sum(fx_edges.*distPerson_edge)/sum(distPerson_edge);
wall_fy = sum(fy_edges.*distPerson_edge)/sum(distPerson_edge);


% Plot of the four points, the person position at the forces applied on the
% person
figure(1)
hold on
plot(person(1),person(2),'ro')
plot(x_edges,y_edges,'b+')
quiver(x_edges,y_edges,fx_edges./(fx_edges.^2+fy_edges.^2).^(1/2),fy_edges./(fx_edges.^2+fy_edges.^2).^(1/2),0.9)
quiver(person(1),person(2),wall_fx/(wall_fx^2+wall_fy^2).^(1/2),wall_fy/(wall_fx^2+wall_fy^2).^(1/2),0.8)
hold off
end



