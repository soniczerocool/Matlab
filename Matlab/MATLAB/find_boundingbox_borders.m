function [ymin,ymax,xmin,xmax,zmin,zmax] = find_boundingbox_borders(T1C)

x_vector = sum(sum(T1C,3),2);
[x_idx x_vector_zero] = find(x_vector>0.01);
xmin = x_idx(1);
xmax = x_idx(end);


y_vector = sum(sum(T1C,3),1);
y_vector = reshape(y_vector,size(T1C,2),1);
[y_idx y_vector_zero] = find(y_vector>0.01);
ymin = y_idx(1);
ymax = y_idx(end);


z_vector = sum(sum(T1C,1),2);
z_vector = reshape(z_vector,size(T1C,3),1);
[z_idx z_vector_zero] = find(z_vector>0.01);
zmin = z_idx(1);
zmax = z_idx(end);
