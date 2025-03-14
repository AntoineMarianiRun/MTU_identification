function plotmodel(muscle_origins, muscle_insertions, joint_centers)
% PLOTMODEL Plots a 3D representation of muscle origins, insertions, and joint centers.
%
% USAGE:
%   plotmodel(muscle_origins, muscle_insertions, joint_centers)
%
% INPUTS:
%   muscle_origins       - A 3xN matrix where each column represents the 3D coordinates of a muscle's origin.
%   muscle_insertions    - A 3xN matrix where each column represents the 3D coordinates of a muscle's insertion.
%   joint_centers        - A 3xN matrix where each column represents the 3D coordinates of a joint center.
%
% DESCRIPTION:
%   This function creates a 3D plot where muscle origins are represented as red dots,
%   muscle insertions as blue dots, and joint centers as black dots. It also draws lines
%   between joint centers and between muscle origins and their respective insertions.
%
% NOTE:
%   It is assumed that the matrices have at least 4 columns. If not, the function might error out.
%
% EXAMPLE:
%   plotmodel(rand(3,4), rand(3,4), rand(3,4));



origins = full(muscle_origins);
insertions = full(muscle_insertions);
markers = full(joint_centers);
%%
figure("Name","model", "Color",[1 1 1])
plot3(origins(1,:), origins(2,:), origins(3,:),'or')
hold on
plot3(insertions(1,:), insertions(2,:), insertions(3,:),'ob')
plot3(markers(1,:), markers(2,:), markers(3,:),'ok')

for i = 1:size(markers,2)-1
    line([markers(1,i),markers(1,i+1)], [markers(2,i),markers(2,i+1)], [markers(3,i),markers(3,i+1)],'Color','black')
end

for i = 1:size(origins,2)
    line([origins(1,i),insertions(1,i)], [origins(2,i),insertions(2,i)], [origins(3,i),insertions(3,i)],'Color','red')
end

legend([{'muscle origin'},{'muscle insersion'},{'joint center'}])
xlim([-1 1])
ylim([-1 1])
axis equal
xlabel('x')
ylabel('y')
zlabel('z')

end
