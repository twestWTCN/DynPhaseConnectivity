function [x y] = remnan(x,y)
x(isnan(y)) = []; y(isnan(y)) = [];
y(isnan(x)) = []; x(isnan(x)) = [];