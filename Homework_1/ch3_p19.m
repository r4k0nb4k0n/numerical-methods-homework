#  2017-2 Mathematical Models for Engineering Problems and Differential Equations
#	Numerical Methods for Engineers and Scientists
#	Chapter 3 Problem 19
#	Written by Choe Hyeong Jin, Dept. of Computer Science, Univ. of Seoul

function Xs = QuadSecRoot(F, a, b)
	imax = 1000; tol = abs((b - a)/2);
	Fa = F(a); Fb = F(b);
	if Fa*Fb > 0
		disp('Error: The function has the same sign at points a and b.');
	else
		for i = 1:imax
			xNS = [0, 0, 0, 0, 0];
			FxNS = [0, 0, 0, 0, 0];
			gap = (b-a)/4;
			for j = 1:5
				xNS(j) = a + gap * (j-1);
				FxNS(j) = F(xNS(j));
			end
			tol = 0.000001*xNS(j);
			toli = (b-a)/2;
			if toli < tol
				break;
			end
			for j = 1:5
				if FxNS(j) == 0
					fprintf('An exact solution x = %11.6f was found\n', xNS(j));
					Xs = xNS(j);
					return;
				end
			end
			for j = 2:5
				if FxNS(j-1) * FxNS(j) < 0
					a = xNS(j-1);
					b = xNS(j);
					break;
				end
			end
			if i == imax
				fprintf('Solution was not obtained in %i iterations\n', imax);
				break;
			end
		end
		Xs = xNS(3);
	end
endfunction
