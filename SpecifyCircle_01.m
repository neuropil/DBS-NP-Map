function [xunit, yunit] = SpecifyCircle_01(x,y,r,numPts)


if nargin == 3
    numPts = 100;
end

steps = round(numPts/2);

th = 0:pi/steps:2*pi;
xunit = r * cos(th) + x;
yunit = r * sin(th) + y;
