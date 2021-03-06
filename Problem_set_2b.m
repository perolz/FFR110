clc, clear all
N=128;
iterations=10000;
u=zeros(N,N);
v=zeros(N,N);
a=3;
b=8;
u0=a;
v0=b/a;
deltat=0.001;
deltax=0.4;

Dv=9;
for i=1:N
    for j=1:N
       r=rand;
       u(i,j)=0.9*u0+0.2*u0*r;
       v(i,j)=0.9*v0+0.2*v0*r;
    end
end
set(gcf, 'Units', 'Normalized', 'OuterPosition', [0, 0.04, 1, 0.96]);
set(gca,'nextplot','replacechildren')
contourf(u,5)
cbr=colorbar;
caxis([2.8 3.3])
xlabel('X1')
ylabel('X2')
k=1;
title(sprintf('u after %d iterations with %c t = %1.3f',k,916,deltat))

%set(cbr,'YTick',2.5:0.5:4)
writerObj = VideoWriter('WaveMovie.avi'); % Name it.
writerObj.FrameRate = 10;
open(writerObj);
frame = getframe(gcf); % 'gcf' can handle if you zoom in to take a movie.
writeVideo(writerObj, frame);

for k=1:iterations
   uOld=u;
   vOld=v;
   for i=1:N
       for j=1:N
           i1=i-1;
           i2=i+1;
           j1=i-1;
           j2=i+1;
           if(i1<1)
               i1=N;
           end
           if(i2>N)
               i2=1;
           end
           if(j1<1)
               j1=N;
           end
           if(j2>N)
               j2=1;
           end
           
           u(i,j)=uOld(i,j)+deltat*(a-(b + 1)*uOld(i,j) + uOld(i,j)^2*vOld(i,j) +...
           (uOld(i2,j)+uOld(i1,j)+uOld(i,j2)+uOld(i,j1)-4*uOld(i,j))/deltax^2);
           v(i,j)=vOld(i,j)+deltat*(b*uOld(i,j)-uOld(i,j)^2*vOld(i,j)+...
               Dv*(vOld(i2,j)+vOld(i1,j)+vOld(i,j2)+vOld(i,j1)-4*vOld(i,j))/deltax^2);
       end

   end
%    if(mod(k,20)==0)
%        contourf(u)
%        title(sprintf('u after %d iterations with %c t = %1.3f',k,916,deltat))
%        frame = getframe(gcf); % 'gcf' can handle if you zoom in to take a movie.
%        writeVideo(writerObj, frame);
%    end
   if(k==1000)
      contourf(u)
      title(sprintf('u after %d iterations with %c t = %1.3f',k,916,deltat))
      saveas(gcf,sprintf('u_2b_%1.1f.png',Dv)) 
   end
end
close(writerObj)
contourf(u)
caxis([2.8 6])
title(sprintf('u after %d iterations with %c t = %1.3f',k,916,deltat))
saveas(gcf,sprintf('u_end_2b_%1.1f.png',Dv))
%%
for j=1:iterations-1
    for i=1:N
        if(i==1)
            u(i,j+1)=a;
            I(i,j+1)=I(i,j)+deltat*(-c*I(i,j)-I(i,j)*(I(i,j)+S(i,j))/K+a*S(i,j)*I(i,j)...
                +D*(I(i+1,j)+I(N,j)-2*I(i,j))/deltax);
        elseif i==N
            S(i,j+1)=S(i,j)+deltat*(b*(I(i,j)+...
                S(i,j))-c*S(i,j)-S(i,j)*(I(i,j)+S(i,j))/K-a*S(i,j)*I(i,j)+D*(S(1,j)+S(i-1,j)-2*S(i,j))/deltax);
            I(i,j+1)=I(i,j)+deltat*(-c*I(i,j)-I(i,j)*(I(i,j)+S(i,j))/K+a*S(i,j)*I(i,j)...
                +D*(I(1,j)+I(i-1,j)-2*I(i,j))/deltax);
        else
            I(i,j+1)=I(i,j)+deltat*(-c*I(i,j)-I(i,j)*(I(i,j)+S(i,j))/K+a*S(i,j)*I(i,j)...
                +D*(I(i+1,j)+I(i-1,j)-2*I(i,j))/deltax);
            S(i,j+1)=S(i,j)+deltat*(b*(I(i,j)+...
                S(i,j))-c*S(i,j)-S(i,j)*(I(i,j)+S(i,j))/K-a*S(i,j)*I(i,j)+D*(S(i+1,j)+S(i-1,j)-2*S(i,j))/deltax);
        end
    end
end


