function [x , y] = get_points(f)


%%% Get Initial Contour
fig = figure;
imshow(f,[]);
[y,x]=getpts(fig);



end