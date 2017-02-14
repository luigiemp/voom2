clear all;
close all;
clc;

% Load Fiber dir
f_vec = load('/u/home/l/luigiemp/project-cardio/CardiacMesh/Contraction/RabbitLV/CoarseLV.fiber');
fileID = fopen('/u/home/l/luigiemp/project-cardio/CardiacMesh/Contraction/RabbitLV/RabbitLV.fiber','w');

for i = 1:size(f_vec,1)
    f = f_vec(i,:);
    f = f/norm(f);
    m(1) = f(2);
    m(2) = -f(1);
    m(3) = 0;
    m = m/norm(m);
    n = cross(f,m);
    
    fprintf(fileID, '%1.8f %1.8f %1.8f %1.8f %1.8f %1.8f %1.8f %1.8f %1.8f\n', ....
            [f,m,n] );
end

fclose(fileID);