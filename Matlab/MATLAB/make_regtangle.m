function [xmin,xmax,ymin,ymax,zmin,zmax] = make_regtangle(T1C)

slice = input('enter slice number to choose the regtangle from:');
fig = T1C(:,:,slice);
imshow(fig,[]);
rect = getrect;
xmin = rect(1);ymin=rect(2); width=rect(3); height=rect(4);
xmax = xmin + width;
ymax = ymin + height;

zmin = input('abnormality start slice number: ');
zmax = input('abnormality end slice number: ');
