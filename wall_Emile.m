clc
close all
clear

%Choose 1 if you want to save the static force field in a file, choose 0 if
%you have already saved it and you just want to access to the saved file
%data.
doYouWantToUsedSavedForces = 1;
%Choose 1 if you want to record a video about the construction of the
%static force field.
%IT TAKES HOURS
doYouWantToMakeAVideo = 0;



% Definition of the (x,y) coordinate system (on bottom left)
x_map = linspace(0,10,200);
y_map = linspace(0,10,200);
% Building the map
[X_map,Y_map] = meshgrid(x_map,y_map(end:-1:1));

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

if doYouWantToUsedSavedForces ~= 0
    [static_fx,static_fy] = readForceFile();
else 
    tic
    for i=1:max(size(wall(1,1,:)))
        [static_fx,static_fy] = staticForceFunction(wall(:,:,i),X_map,Y_map,static_fx,static_fy,distanceAmongPointsOnTheWall,wallRepulsionConstant);
    end
    toc
    saveForcesInFiles(static_fx,static_fy);
    
end



if doYouWantToMakeAVideo == 1
    static_fx = zeros(max(size(x_map)));
    static_fy = zeros(max(size(x_map)));
    writerObj = VideoWriter('video.avi');
    writerObj.FrameRate = 60;
    open(writerObj);
    for i=1:max(size(wall(1,1,:)))
        [static_fx,static_fy] = recordVideo(writerObj,wall(:,:,i),X_map,Y_map,static_fx,static_fy,distanceAmongPointsOnTheWall,wallRepulsionConstant);

    end
    close(writerObj);
end


coloringTheMap(static_fx,static_fy,distanceAmongPointsOnTheWall,X_map,Y_map,wallRepulsionConstant)
for i=1:max(size(wall(1,1,:)))
    wallPlotter(wall(:,:,i));
end
functionToPlotTheStaticField(X_map,Y_map,static_fx,static_fy)

% Example of a person in the map
person = [3.2 7.5];
% Calculation (and plot) of the forces applied by wall on that person
tic
[wall_fx,wall_fy]=person_wallInteraction(person,static_fx,static_fy,X_map,Y_map);
toc














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
    static_fy = static_fy + wallRepulsionConstant.*(Y_map-y_wall(k))./((X_map-x_wall(k)).^2+(Y_map-y_wall(k)).^2).^(3/2);
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

%{
% Plot of the four points, the person position at the forces applied on the
% person
figure(1)
hold on
plot(person(1),person(2),'bo')
plot(x_edges,y_edges,'r+')
q = quiver(x_edges,y_edges,fx_edges./(fx_edges.^2+fy_edges.^2).^(1/2),fy_edges./(fx_edges.^2+fy_edges.^2).^(1/2),3.5);
q.Color = 'red';
q.LineWidth  =0.9;
q = quiver(person(1),person(2),wall_fx/(wall_fx^2+wall_fy^2).^(1/2),wall_fy/(wall_fx^2+wall_fy^2).^(1/2),0.3);
q.Color = 'yellow';
q.LineWidth  = 2;
q.MaxHeadSize = 0.6;

hold off
%}
end

function saveForcesInFiles(static_fx,static_fy)
fileID = fopen('wall_fx.txt','w');
fprintf(fileID,'%f\n',static_fx);
fclose(fileID);

fileID = fopen('wall_fy.txt','w');
fprintf(fileID,'%f\n',static_fy);
fclose(fileID);
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
heatMap = (statif_tot_f./Fmax).^(1/3);


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

function wallPlotter(wall)
%Plot of the discretized wall
figure(1)
hold on
plot(wall(:,1),wall(:,2),'r')
hold off
end

function [static_fx,static_fy] = recordVideo(writerObj,wall,X_map,Y_map,static_fx,static_fy,distanceAmongPointsOnTheWall,wallRepulsionConstant)

%SAME CODE FOR THE CALCULATION OF THE FORCE
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
wall_fx = zeros(size(static_fx));
wall_fy = zeros(size(static_fy));
for k=1:numberOfPointsOnTheWall
    wall_fx = wall_fx + wallRepulsionConstant.*(X_map-x_wall(k))./((X_map-x_wall(k)).^2+(Y_map-y_wall(k)).^2).^(3/2);
    wall_fy = wall_fy + wallRepulsionConstant.*(Y_map-y_wall(k))./((X_map-x_wall(k)).^2+(Y_map-y_wall(k)).^2).^(3/2);
end






line =@(y0,x0,m,x) y0 + m.*(x-x0);

%Initialization of the map tracking the distance of every point from the
%wall
distWall = zeros(size(X_map));
%m must be different from infinite
if wall(1,2)~=wall(2,2)
    y1 = max(wall(:,2));
    x1 = (wall(:,2)==y1)'*[wall(1,1);wall(2,1)];
    y0 = min(wall(:,2));
    x0 = (wall(:,2)==y0)'*[wall(1,1);wall(2,1)];
    
    if wall(1,1)~=wall(2,1)
        % Angular coefficient of the perpendicular line to the wall
        m = (wall(1,2)-wall(2,2))/(wall(1,1)-wall(2,1));
        mp = - 1/m;
        distWall(distWall==0) = abs(wall(1,2)-m*wall(1,1)+m.*X_map(distWall==0)-Y_map(distWall==0))/sqrt(1+m^2);
    else
        mp = 0;
        distWall(Y_map<=y1 & Y_map>=y0) = abs(X_map(Y_map<=y1 & Y_map>=y0)-x0);
    end
    distWall(Y_map>line(y1,x1,mp,X_map)) = sqrt((X_map(Y_map>line(y1,x1,mp,X_map)) -x1).^2+(Y_map(Y_map>line(y1,x1,mp,X_map)) -y1).^2);
    distWall(Y_map<line(y0,x0,mp,X_map)) = sqrt((X_map(Y_map<line(y0,x0,mp,X_map)) -x0).^2+(Y_map(Y_map<line(y0,x0,mp,X_map)) -y0).^2);
    %surf(X_map,Y_map,distWall)
else
    x1 = max(wall(:,1));
    y1 = (wall(:,1)==x1)'*[wall(1,2);wall(2,2)];
    x0 = min(wall(:,1));
    y0 = (wall(:,1)==x1)'*[wall(1,2);wall(2,2)];
    distWall(X_map<x0) = sqrt((X_map(X_map<x0) -x0).^2+(Y_map(X_map<x0) -y0).^2);
    distWall(X_map>x1) = sqrt((X_map(X_map>x1) -x1).^2+(Y_map(X_map>x1) -y1).^2);
    distWall(distWall==0) = abs(Y_map(distWall==0) -y0);
    %surf(X_map,Y_map,distWall)

end


distMin = distanceAmongPointsOnTheWall;
distMax = max(max(X_map));

distWall(distWall>distMax) = distMax;

updatedFx = static_fx;
updatedFy = static_fy;
fx = updatedFx./sqrt(updatedFx.^2+updatedFy.^2);
fy = updatedFy./sqrt(updatedFx.^2+updatedFy.^2);
figure(1);
coloringTheMap(updatedFx,updatedFy,distanceAmongPointsOnTheWall,X_map,Y_map,wallRepulsionConstant)
wallPlotter(wall);
hold on
q = quiver(X_map(1:3:end,1:3:end),Y_map(1:3:end,1:3:end),fx(1:3:end,1:3:end),fy(1:3:end,1:3:end),0.6);
q.Color = 'red';
hold off
axis([0 max(max(X_map)) 0 max(max(Y_map))]);
frame = getframe(gcf);
writeVideo(writerObj,frame);
pause(0.7)


for i=10:-0.05:1
    if((distMax).^(0.6)/i>=max(max(distWall)))
        continue
    end
    updatedFx = static_fx;
    updatedFy = static_fy;
    updatedFx((distWall).^(0.6)<=(distMax).^(0.6)/i) = static_fx((distWall).^(0.6)<=(distMax).^(0.6)/i)+wall_fx((distWall).^(0.6)<=(distMax).^(0.6)/i);
    updatedFy((distWall).^(0.6)<=(distMax).^(0.6)/i) = static_fy((distWall).^(0.6)<=(distMax).^(0.6)/i)+wall_fy((distWall).^(0.6)<=(distMax).^(0.6)/i);
    fx = updatedFx./sqrt(updatedFx.^2+updatedFy.^2);
    fy = updatedFy./sqrt(updatedFx.^2+updatedFy.^2);
    
    figure(1);
    coloringTheMap(updatedFx,updatedFy,distanceAmongPointsOnTheWall,X_map,Y_map,wallRepulsionConstant)
    wallPlotter(wall);
    hold on
    q = quiver(X_map(1:3:end,1:3:end),Y_map(1:3:end,1:3:end),fx(1:3:end,1:3:end),fy(1:3:end,1:3:end),0.6);
    q.Color = 'red';
    hold off
    axis([0 max(max(X_map)) 0 max(max(Y_map))]);
    frame = getframe(gcf);
    writeVideo(writerObj,frame);
    pause(0.1)
end
static_fx = static_fx + wall_fx;
static_fy = static_fy + wall_fy;
end
