#  2017-2 Mathematical Models for Engineering Problems and Differential Equations
#	Numerical Methods for Engineers and Scientists
#	Chapter 4 Problem 27
#	Written by Choe Hyeong Jin, Dept. of Computer Science, Univ. of Seoul

function c = CondNumb_One(A)
	c = norm(A, 1) * norm(inv(A), 1);
endfunction