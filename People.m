% this file is a sample program to illustract the structure of our social
% model algorithm for Apero

N_p = 20;   % number of people
N_t = 10;  % number of position
N_f = 2;   % number of foods

position_people = rand(N_p, 2);  % the x (first column), y (second column) coordinates of position of people (initialization)
position_table = rand(N_t, 2);  % the position of tables
position_food = rand(N_f, 2);  % position of foods

function [d] = distance(people_1, people_2)
d = sqrt((people_1(1)-people_2(1)).^2 + (people_1(2)-people_2(2)).^2);
end

function [force] = force(people_1, people_2)

end

function [force] = force(people, table)
end

function [force] = force(people, food)
end

function [force] = force(people, wall)
end


