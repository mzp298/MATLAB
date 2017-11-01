clear all

clc

t = -2*pi:0.1:2*pi;

figure('name','hello')

filename = 'sin2x.gif';

n=length(t);

y = sin(2*t);

plot(t,y,'b')

hold on

for i = 1:n

   xx = t(i);

   yy = y(i);

   plot(xx,yy,'r.','markersize',16);

   title('sin(2x)')

   xlabel('x')

   ylabel('y')

   axis([-2*pi,2*pi,-2,2])

   grid on

   frame = getframe(1);

   im = frame2im(frame);

   [A,map] = rgb2ind(im,256);

   if i == 1;

       imwrite(A,map,filename,'gif','LoopCount',Inf,'DelayTime',0.01);

   else

       imwrite(A,map,filename,'gif','WriteMode','append','DelayTime',0.01);

   end

end