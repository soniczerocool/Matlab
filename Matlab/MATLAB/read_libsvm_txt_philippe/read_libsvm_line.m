function [label x y z] = read_libsvm_line(line)
% takes in a line as libsvm file format as 'label 1:x 2:y 3:z' and gives out [label x y z
[label x_str y_str z_str] = strread(line, '%s %s %s %s', 'delimiter', ' ');
[ix x] = strtok(x_str,'0');
[iy y] = strtok(y_str, '0');
[iz z] = strtok(z_str, '0');
x = str2double(x);
y = str2double(y);
z = str2double(z);
label = str2double(label);
