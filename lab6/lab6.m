%% Laplace Transform Lab: Solving ODEs using Laplace Transform in MATLAB
% This lab will teach you to solve ODEs using a built in MATLAB Laplace transform 
% function |laplace|.
% 
% There are five (5) exercises in this lab that are to be handed in. Write your 
% solutions in a separate file, including appropriate descriptions in each step.
% 
% Include your name and student number in the submitted file.

%% Student Information
% 
%  Student Name: Linda Zhao
%
%  Student Number: 1008107683
%

%% Using symbolic variables to define functions
% In this exercise we will use symbolic variables and functions.

syms t s x y

f = cos(t)
h = exp(2*x)

%% Laplace transform and its inverse

% The routine |laplace| computes the Laplace transform of a function

F=laplace(f)

% By default it uses the variable |s| for the Laplace transform But we can specify 
% which variable we want:

H=laplace(h)
laplace(h,y)

% Observe that the results are identical: one in the variable |s| and the
% other in the variable |y|

% We can also specify which variable to use to compute the Laplace transform:

j = exp(x*t)
laplace(j)
laplace(j,x,s)

% By default, MATLAB assumes that the Laplace transform is to be computed
% using the variable |t|, unless we specify that we should use the variable
% |x|

% We can also use inline functions with |laplace|. When using inline functions, 
% we always have to specify the variable of the function.

l = @(t) t^2+t+1
laplace(l(t))

% MATLAB also has the routine |ilaplace| to compute the inverse Laplace transform

ilaplace(F)
ilaplace(H)
ilaplace(laplace(f))

% If |laplace| cannot compute the Laplace transform, it returns an unevaluated 
% call.

g = 1/sqrt(t^2+1)
G = laplace(g)

%% % But MATLAB "knows" that it is supposed to be a Laplace transform of a function. 
% So if we compute the inverse Laplace transform, we obtain the original function

ilaplace(G)

% The Laplace transform of a function is related to the Laplace transform of 
% its derivative:

syms g(t)
laplace(diff(g,t),t,s)

%% Exercise 1
% Objective: Compute the Laplace transform and use it to show that MATLAB 'knows' 
% some of its properties.
% 
% Details: 
% 
% (a) Define the function |f(t)=exp(2t)*t^3|, and compute its Laplace transform 
% |F(s)|. (b) Find a function |f(t)| such that its Laplace transform is  |(s - 
% 1)*(s - 2))/(s*(s + 2)*(s - 3)| (c) Show that MATLAB 'knows' that if |F(s)| 
% is the Laplace transform of  |f(t)|, then the Laplace transform of |exp(at)f(t)| 
% is |F(s-a)| 
% 
% (in your answer, explain part (c) using comments). 
% 
% Observe that MATLAB splits the rational function automatically when solving 
% the inverse Laplace transform.

% ============================================================================
% Exercise 1 Submission
% ============================================================================

% Define symbolic variables
syms t s

% (a)
% Define inline function
f(t) = exp(2*t)*(t^3);
% Compute (and display) laplace transform
F = laplace(f)

% (b)
% Define Laplace transformed function
G(s) = ((s-1)*(s-2))/(s*(s+2)*(s-3));
% Compute (and display) inverse
g = ilaplace(G)

% (c)
syms q(t) p(t) a t

% Define function
q(t) = exp(a*t)*p(t);
% Compute laplace transforms
P = laplace(p)
Q = laplace(q)

% When displaying the laplace transforms of the arbitrary functions p(t) and
% q(t) = exp(a*t)*p(t) (represented by P and Q, respectively), 
% MATLAB simply changes the last argument of the |laplace| routine from 
% |s| to |s-a|. In other words, MATLAB simply replaces the |s| parameter in
% P(s) with |s-a| to evaluate Q, meaning that MATLAB "knows" we can write 
% Q = P(s-a).

% ============================================================================

%% Heaviside and Dirac functions
% These two functions are builtin to MATLAB: |heaviside| is the Heaviside function 
% |u_0(t)| at |0|
% 
% To define |u_2(t)|, we need to write

f=heaviside(t-2)
ezplot(f,[-1,5])

% The Dirac delta function (at |0|) is also defined with the routine |dirac|

g = dirac(t-3)

% MATLAB "knows" how to compute the Laplace transform of these functions

laplace(f)
laplace(g)

%% Exercise 2
% Objective: Find a formula comparing the Laplace transform of a translation 
% of |f(t)| by |t-a| with the Laplace transform of |f(t)|
% 
% Details: 

% * Give a value to |a|
% * Let |G(s)| be the Laplace transform of |g(t)=u_a(t)f(t-a)| and |F(s)| is 
% the Laplace transform of |f(t)|, then find a formula relating |G(s)| and |F(s)|

% In your answer, explain the 'proof' using comments.

% ============================================================================
% Exercise 2 Submission
% ============================================================================

syms u_a(t) f(t) g(t) t

for i = 1:5
    a = i
    % Define Heaviside function u_a(t) 
    u_a(t) = heaviside(t-a);

    % Define function g(t)
    g(t) = u_a(t)*f(t-a);

    % Compute and display Laplace transforms for several values of a
    F = laplace(f);
    G = laplace(g)
end

% As shown by the MATLAB output display, the Laplace transform G of g(t) =
% u_a(t)*f(t-a) is simply exp(-a*s)*F, where F is the Laplace transform of
% f(t) and a E [1,5] in this example.
% This outcome is corroborated by Theorem 5.5.1 in the textbook,
% which states that L{u_a(t)*f(t-a)} = exp(-a*s)*L{f(t)}. The proof of this
% theorem lies in a change of variables used to transform the integral
% (from c to infinity) of exp(-st)*f(t-c)dt to the integral of
% exp(-(z+c)*s)*f(x)dz, where z = t - c. Then, the exp(-c*s) term can be
% taken out of the integral since the integral is with respect to z, which
% results in L{u_a(t)*f(t-a)} = exp(-c*s) * integral from 0 to infinity of
% exp(-sz)*f(z)dz = exp(-c*s)*F(s). 
% ============================================================================

%% Solving IVPs using Laplace transforms
% Consider the following IVP, |y''-3y = 5t| with the initial conditions |y(0)=1| 
% and |y'(0)=2|. We can use MATLAB to solve this problem using Laplace transforms:

% First we define the unknown function and its variable and the Laplace
% tranform of the unknown

syms y(t) t Y s

% Then we define the ODE

ODE=diff(y(t),t,2)-3*y(t)-5*t == 0

% Now we compute the Laplace transform of the ODE.

L_ODE = laplace(ODE)

% Use the initial conditions

L_ODE=subs(L_ODE,y(0),1)
L_ODE=subs(L_ODE,subs(diff(y(t), t), t, 0),2)

% We then need to factor out the Laplace transform of |y(t)|

L_ODE = subs(L_ODE,laplace(y(t), t, s), Y)
Y=solve(L_ODE,Y)

% We now need to use the inverse Laplace transform to obtain the solution
% to the original IVP

y = ilaplace(Y)

% We can plot the solution

ezplot(y,[0,20])

% We can check that this is indeed the solution

diff(y,t,2)-3*y
%% Exercise 3
% Objective: Solve an IVP using the Laplace transform
% 
% Details: Explain your steps using comments

% * Solve the IVP
% * |y'''+2y''+y'+2*y=-cos(t)|
% * |y(0)=0|, |y'(0)=0|, and |y''(0)=0|
% * for |t| in |[0,10*pi]|
% * Is there an initial condition for which |y| remains bounded as |t| goes 
% to infinity? If so, find it.

% ============================================================================
% Exercise 3 Submission
% ============================================================================

% Define unknown function, independent variable and Laplace transform of
% the unknown function 
syms y(t) t Y s

% Define ODE
ODE = diff(y(t),t,3) + 2*diff(y(t),t,2) + diff(y(t),t,1) + 2*y(t) + cos(t) == 0;

% Compute Laplace transform of ODE
L_ODE = laplace(ODE);

% Use initial conditions
L_ODE = subs(L_ODE,y(0),0);
L_ODE = subs(L_ODE,subs(diff(y(t),t),t,0),0);
L_ODE = subs(L_ODE,subs(diff(y(t),t,2),t,0),0);

% Factor out Laplace transform of y(t) and solve for it
L_ODE = subs(L_ODE,laplace(y(t), t, s), Y)
Y = solve(L_ODE,Y)

% Get (and display) the solution to the original IVP using the inverse Laplace transform 
y = ilaplace(Y)

% Plot solution
ezplot(y,[0,10*pi]);
title('Solution to Exercise 3 IVP')
xlabel('t');
ylabel('y');

% Based on the solution to the problem, there is no intial condition for
% which y remains bounded as t goes to infinity. This is because the
% solution contains terms containing t multiplied by a sinusoidal function,
% which will always oscillate while increasing to infinity as t goes to
% infinity. Changing the initial conditions will not change this fundamental 
% form of the solution.
% ============================================================================

%% Exercise 4
% Objective: Solve an IVP using the Laplace transform
% 
% Details: 

% * Define 
% * |g(t) = 3 if 0 < t < 2|
% * |g(t) = t+1 if 2 < t < 5|
% * |g(t) = 5 if t > 5|
% * Solve the IVP
% * |y''+2y'+5y=g(t)|
% * |y(0)=2 and y'(0)=1|
% * Plot the solution for |t| in |[0,12]| and |y| in |[0,2.25]|.

% In your answer, explain your steps using comments.

% ============================================================================
% Exercise 4 Submission
% ============================================================================

% Define unknown function, independent variable and Laplace transform of
% the unknown function 
syms y(t) t Y s

% Define heaviside functions for g(t)
u_0(t) = heaviside(t);
u_2(t) = heaviside(t-2);
u_5(t) = heaviside(t-5);

% Define g(t) using heaviside functions, whose "coefficients" can be found 
% by graphing each piecewise component.
g(t) = 3*u_0(t) + (t-2)*u_2(t) + (4-t)*u_5(t);

% Define ODE
ODE = diff(y(t),t,2) + 2*diff(y(t),t,1) + 5*y(t) - g(t) == 0;

% Compute Laplace transform of ODE
L_ODE = laplace(ODE);

% Use initial conditions
L_ODE = subs(L_ODE,y(0),2);
L_ODE = subs(L_ODE,subs(diff(y(t),t),t,0),1);

% Factor out Laplace transform of y(t) and solve for it
L_ODE = subs(L_ODE,laplace(y(t), t, s), Y)
Y = solve(L_ODE,Y)

% Get (and display) the solution to the original IVP using the inverse Laplace transform 
y = ilaplace(Y)

% Plot solution
ezplot(y,[0,12,0,2.25]);
title('Solution to Exercise 4 IVP')
xlabel('t');
ylabel('y');

% ============================================================================

%% Exercise 5
% Objective: Use the Laplace transform to solve an integral equation
% 
% Verify that MATLAB knowns about the convolution theorem by explaining why 
% the following transform is computed correctly.

% ============================================================================
% Exercise 5 Submission
% ============================================================================

syms t tau y(tau) s
I=int(exp(-2*(t-tau))*y(tau),tau,0,t)
laplace(I,t,s)

% The Convolution Theorem (5.8.3) states that for two functions f and g, 
% the Laplace transform of f convolved with g is equal to the product of
% the individual Laplace transforms of f and g (i.e. L{h(t)} = F(s)*G(s)
% where h(t) = integral from 0 to t of f(t-tau)*g(tau)d(tau)).

% In the case of the above transform, we can take f(t) = exp(-2*t) as convolved
% with some arbitrary function g(t) = y(t). The laplace transform of this
% convolution is taken in the following line and the result is displayed.

% We can prove that MATLAB computes the transform correctly by writing out
% explicitly the product of the two individual Laplace transforms:

syms f(t) g(t)

f(t) = exp(-2*t);

F = laplace(f)
G = laplace(y)

H = F * G

% Since H(s) = F(s)*G(s) is the same as the above output (given by laplace(I,t,s), 
% by the convolution theorem, we know that MATLAB computed the transform correctly.

% ============================================================================