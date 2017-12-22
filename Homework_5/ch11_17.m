% ch11_17.m
% Written by Jeong Ye Ji
% Commented by Choe Hyeong Jin

clear all;

function [x,y] = BVP2ndDriv(pOFx, qOFx, rOFx, a, b, Ya, Db, n)
  h=(b-a)/n;
  x=a:h:b;
  y=zeros(1,n+1);
  ys=zeros(1,n-1);
  c=zeros(1,n-1);
  a=zeros(n-1,n-1);

  for i=1:n-1
    c(1,i)=rOFx(x(i))*2*h*h;
  end

  a(1,1)=-(4-2*h*h*qOFx(x(2)));
  a(1,2)=2+h*pOFx(x(2));
  c(1,1)=c(1,1)-(2-h*pOFx(x(2)));
    
  for i=2:n-2
    a(i,i-1)=2-h*pOFx(x(i+1));
    a(i,i)=-(4-2*h*h*qOFx(x(i+1)));
    a(i,i+1)=2+h*pOFx(x(i+1));
  end
    
  a(n-1,n-2) = 2-h*pOFx(x(n));
  a(n-1,n-1) = -(4-2*h*h*qOFx(x(n)));
  c(1,n-1) = c(1,n-1)-(2+h*pOFx(x(n)));

  ys=c*inv(a);

  y(1,1)=Ya;
  for i=2:n
    y(1,i)=ys(1,i-1);
  end
  y(1,n+1) = 2*h*Db + y(1,49);    
end

[x,y] = BVP2ndDriv(@(x) 1/x,@(x) 0,@(x) -10, 1, 3, 1, -1.2, 50);
plot(x,y);
title("BVP2ndDriv()");
xlabel("x");
ylabel("y");
