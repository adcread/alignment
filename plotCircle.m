function [  ] = plotCircle( Re, Im, r )
%PLOTCIRCLE Summary of this function goes here

ang=0:0.01:2*pi; 
xp=r*cos(ang);
yp=r*sin(ang);
plot(Re+xp,(Im+yp));

end

