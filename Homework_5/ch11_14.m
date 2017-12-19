function [x, t, y] = Sys2ODEsRK4(ODE1, ODE2, a, b, h, x1, y1)
% This solves a system of two first-order inital value ODEs using
% fourth-order Runge-Kutta method.
% The independent variable is t, and the dependent variables are x, y.
% Input variables : 
% ODE1 -> dx/dt
% ODE2 -> dy/dt
% a    -> The first value of t.
% b    -> The last value of t.
% h    -> The size of a increment.
% x1   -> The initial value of x.
% y1   -> The initial value of y.
% Output variables : 
% x    -> A vector with the x coordinates of the solution points.
% t    -> A vector with the t coordinates of the solution points.
% y    -> A vector with the y coordinates of the solution points.
  t(1) = a; x(1) = x1; y(1) = y1;
  n = (b - a)/h;
  for i = 1:n
    t(i+1) = t(i) + h;
    tm = t(i) + h/2;
    Kx1 = ODE1(t(i), x(i), y(i));
    Ky1 = ODE2(t(i), x(i), y(i));
    Kx2 = ODE1(tm, x(i) + Kx1*h/2, y(i) + Ky1*h/2);
    Ky2 = ODE2(tm, x(i) + Kx1*h/2, y(i) + Ky1*h/2);
    Kx3 = ODE1(tm, x(i) + Kx2*h/2, y(i) + Ky2*h/2);
    Ky3 = ODE2(tm, x(i) + Kx2*h/2, y(i) + Ky2*h/2);
    Kx4 = ODE1(t(i + 1), x(i) + Kx3*h, y(i) + Ky3*h);
    Ky4 = ODE2(t(i + 1), x(i) + Kx3*h, y(i) + Ky3*h);
    x(i+1) = x(i) + (Kx1 + 2*Kx2 + 2*Kx3 + Kx4)*h/6;
    y(i+1) = y(i) + (Ky1 + 2*Ky2 + 2*Ky3 + Ky4)*h/6;
  end
end

function dwdx = DwDx(x, w, y)
  dwdx = -2*x*w - 5*y + cos(3*x);
end

function dydx = DyDx(x, w, y)
  dydx = w;
end

function [x, y] = BVPShootSecant(fOFx, gOFx, hOfx, a, b, n, Ya, Yb, WL, WH)
% This solves a second-order boundary value problem of the specific form.
% Input variables : 
% fOFx, gOFx, hOFx -> f(x), g(x), h(x)
% a, b             -> defines the domain of the solution
% n                -> The number of subintervals
% Ya, Yb           -> The boundary conditions
% WL, WH           -> The assumed slopes at x = a that are used in the first two solutions
% Output variables :
% x, y             -> The solution of the second-order BVP that is given.

  tol = 0.001;
  h = (b - a)/n;
  [x0, w0, y0] = Sys2ODEsRK4(@DwDx, @DyDx, a, b, h, 0, Ya);
  [x1, w1, y1] = Sys2ODEsRK4(@DwDx, @DyDx, a, b, h, 0, Yb);
  for i = 1:n % Do n times.
    xi = x1 - (x1 - x0).*y1/(y1 - y0) % secant formula
    yi = y1 - (y1 - y0).*x1/(x1 - x0) % secant formula
    if abs((xi - x1)/x1)  < tol
      x = xi;
      y = yi;
      break
    end
    [x0, w0, y0] = Sys2ODEsRK4(@DwDx, @DyDx, a, b, h, xi(n), yi(n));
  end
  if i>n
    fprintf('Solution was not obtained in %i iterations.\n', n);
  else
    plot(x, y)
    xlabel('x'); ylabel('y');
  end
end

[x, y] = BVPShootSecant(@(x) 2*x, 5, @(x) cos(3*x), 0, pi, 10, 1.5, 0, -5, -1.5);
display(x);
