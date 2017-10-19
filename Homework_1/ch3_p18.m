#  2017-2 Mathematical Models for Engineering Problems and Differential Equations
#	Numerical Methods for Engineers and Scientists
#	Chapter 3 Problem 18
#	Written by Choe Hyeong Jin, Dept. of Computer Science, Univ. of Seoul

function X = Ln (p)
	if p <= 0
		disp('Error: p is zero or a negative number.');
		return;
	end
	F = @ (x) exp(x) - p;
	if p > exp(1)
		a = exp(0); b = p;
	else
		a = -1/p; b = exp(0);
	end
	imax = 100; tol = 0.000001;
	Fa = F(a); Fb = F(b);
	if Fa*Fb > 0
		disp('Error: The function has the same sign at points a and b.');
	else
# disp('iteration        a           b    (xNS) Solution    f(xNS)    Tolerance');
		for i = 1:imax
			xNS = (a + b)/2;
			toli = (b - a)/2;
			FxNS = F(xNS);
# fprintf('%9i %11.6f %11.6f %11.6f %11.6f %11.6f\n', i, a, b, xNS, FxNS, toli);
			if FxNS == 0
#	fprintf('An exact solution x = %11.6f was found\n', xNS);
				break;
			end
			if toli < tol
				break;
			end
			if i > imax
				fprintf('Solution was not obtained in %i iterations\n', imax);
				break;
			end
			if F(a)*FxNS < 0
				b = xNS;
			else
				a = xNS;
			end
		end
	end
	X = xNS;
#	fprintf('X = %11.6f\n',X);
endfunction
