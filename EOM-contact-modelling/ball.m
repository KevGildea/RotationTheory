function xdot  = ball(t,x)

global g k c m h r

xdot  = zeros(2,1);

%determine if penetration has occured

pen = (x(2) + r) - h; %compute penetetration for contact

if pen <0 
    xdot(1) = g; xdot(2)  = x(1);
else xdot(1) = g - ((k*pen)/m) -((c*x(1))/m); xdot(2)  = x(1);

end
