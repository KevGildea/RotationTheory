close all; clear

global g k c m h r

g = 9.81;k = 500; c=0; m = 0.2; h=0.2; r = 0.01; 

% set integration time interval, and maximum time step
tspan = [0 2];
tstep = 1e-3;

% set initial conditions: velocity and displacement
X0 = [0  0];

% call ode23 function and feed relevant inputs
OPTIONS = odeset('MaxStep',tstep); %set max integration timestep
[tout,xout] = ode45('ball',tspan,X0,OPTIONS);

% plot position and velocity time history
figure(1); subplot(2,1,1)
set ( gca, 'ydir', 'reverse' )
hold on; grid on; 
plot(tout,xout(:,2),'k.-')
xlabel('time -s')
ylabel('disp - m')
title('ball displacement')
  
figure(1); subplot(2,1,2)
set ( gca, 'ydir', 'reverse' )
hold on; grid on; 
plot(tout,xout(:,1),'k.-')
xlabel('time -s')
ylabel('vel - m/s')
title('ball velocity')

saveas(1,'k5000c0.png')

% Create GIF of ball bouncing
figure(2);
filename = 'k500c0.gif';
x=xout(:,2);
t = tspan(1):tstep:tspan(2);
for n = 1:10:length(x)
      txt = sprintf('t= %d seconds', t(n));
      text(r*2,0,txt) 
      viscircles([0,x(n)],[r],'Color','k')
      xlim([-0.05 0.4]);
      ylim([-0.05 0.4]);
      yline(h, '-','Ground plane', 'LineWidth', 2);
      set ( gca, 'ydir', 'reverse' )
      axis equal
      frame = getframe(2);
      clf
      im = frame2im(frame);
      [imind,cm] = rgb2ind(im,256);
      if n == 1;
          imwrite(imind,cm,filename,'gif', 'DelayTime', 0, 'Loopcount',inf);
      else
          imwrite(imind,cm,filename,'gif', 'DelayTime', 0,'WriteMode','append');
      end
end

