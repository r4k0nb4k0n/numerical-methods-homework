% ch11_14.m
% Written by Lee Tae Hee.
% Commented by Choe Hyeong Jin.

clear
function [x,y,w]=Sys2ODEsRK4(ODE1,ODE2,a,b,h,y1,w1)
  x(1)=a;y(1)=y1;w(1)=w1;
  n=(b-a)/h;  
  for i=1:n
    x(i+1)=x(i)+h;
    xm=x(i)+h/2;
    Ky1=ODE1(x(i),y(i),w(i));
    Kw1=ODE2(x(i),y(i),w(i));
    Ky2=ODE1(xm,y(i)+Ky1*h/2,w(i)+Kw1*h/2);
    Kw2=ODE2(xm,y(i)+Ky1*h/2,w(i)+Kw1*h/2);
    Ky3=ODE1(xm,y(i)+Ky2*h/2,w(i)+Kw2*h/2);
    Kw3=ODE2(xm,y(i)+Ky2*h/2,w(i)+Kw2*h/2);
    Ky4=ODE1(x(i+1),y(i)+Ky3*h,w(i)+Kw3*h);
    Kw4=ODE2(x(i+1),y(i)+Ky3*h,w(i)+Kw3*h);
    y(i+1)=y(i)+(Ky1+2*Ky2+2*Ky3+Ky4)*h/6;
    w(i+1)=w(i)+(Kw1+2*Kw2+2*Kw3+Kw4)*h/6;
  end
end

function [x,y]=BVPShootSecant(fOFx,gOFx,hOFx,a,b,n,Ya,Yb,WL,WH)
  tol = 0.001

  % We have to transform BVP to 2 IVPs.
  ODE1=@(x,y,w) w;
  ODE2=@(x,y,w) hOFx(x)-fOFx(x)*w-gOFx(x)*y; 
   
  % Set the values of w.
  w(1)=WL;
  w(2)=WH;
  
    
  for i=1:n
    [x1,y1,w1]=Sys2ODEsRK4(ODE1,ODE2,a,b,(b-a)/n,Ya,w(i)); % Solve the IVP using Sys2ODEsRK4().
    if(abs(y1(n+1)-Yb)<tol) % Handle the error handling part.
      break;
    end  
    [x2,y2,w2]=Sys2ODEsRK4(ODE1,ODE2,a,b,(b-a)/n,Ya,w(i+1)); % Solve the IVP using Sys2ODEsRK4().
    w(i+2)=w(i+1)-((w(i+1)-w(i))*y2(n+1))/(y2(n+1)-y1(n+1)); % Using the Secant method.
  end
  % We got the solution. 
  x=x1; y=y1;
  plot(x,y);
  title("BVPShootSecant()");
  xlabel("x");
  ylabel("y");
end
[x,y]=BVPShootSecant(@(x) 2*x,@(x) 5, @(x) cos(3*x),0,3,100,1.5,0,-5,-1.5); # pi -> 3
