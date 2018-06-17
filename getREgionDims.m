function [ eQdiameter ] = getREgionDims(x_vec, y_vec)
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here


minX = min(x_vec);
maxX = max(x_vec);

Xdist = maxX - minX;

minY = min(y_vec);
maxY = max(y_vec);

Ydist = maxY - minY;

eQdiameter = mean([Xdist , Ydist]);


end

