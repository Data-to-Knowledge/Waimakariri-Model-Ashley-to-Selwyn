# -*- coding: utf-8 -*-
"""
Author: matth
Date Created: 19/03/2018 12:43 PM
"""

from __future__ import division
from core import env
"""
below if from hydraulics package of stream depletion v3.xls (ecan)
'This module contains routines to calculate groundwater flow solutions.

'This calculates drawdowns for flow to a well from the Theis solution.
'All variables are dimensionless and are defined as follows:
'        W_1=Abs(h)*T/Q     r' = r / L
t
'=tT/(SL^2)
'where Abs(h)=drawdown, Q=well flow, r=radial distance from the pumped well,
't=time, T=transmissivity, S=storage coefficient and L=length that may be
'chosen arbitrarily since it cancels out in the Theis solution. For example,
'choose L=one metre if units of metre are used in the problem. All input and
'output variables are dimensionless, and the prime superscript has been omitted
'for notational convenience.
Function
W_1(r, t)
If
t <= 0
Then
W_1 = 0  #
Else
Pi = 3.141592654
W_1 = Exp1(r ^ 2 / (4 * t)) / (4 * Pi)
End
If
End
Function

'This calculates drawdowns for flow to a well in a leaky aquifer from
'the Hantush solution. All variables are dimensionless and are defined
'as follows:
'    W_2=Abs(h)*T/Q    r' = r / L
t
'=tT/(SL^2)   K=(K' / B')L^2/T
'where Abs(h)=drawdown, Q=well flow, r=radial distance from the pumped well,
't=time, T=transmissivity, S=storage coefficient, K' = aquitard
permeability,
'B' = aquitard
thickness and L = length
that
may
be
chosen
arbitrarily
since
it
'cancels out in the Hantush solution. For example, choose L=one metre if units
'of metre are used in the problem. All input and output variables are dimensionless,
'and the prime superscript has been omitted for notational convenience.
Function
W_2(r, t, K)
If
t <= 0
Then
W_2 = 0  #
Else
Pi = 3.141592654
W_2 = w(r ^ 2 / (4 * t), r * Sqr(K)) / (4 * Pi)
End
IfEnd
Function

'This calculates drawdowns for flow to a well with delayed yield, as first
'conceived by Boulton. All input and output variables are dimensionless and are
'defined As follows:
'   W_3=Abs(h)*T/Q    r' = r / L
t
'=t*T/SL^2  epsilon=S/sigma   K' = (K
"/B")L ^ 2 / T
'where Abs(h)=drawdown, Q=well flow, T=transmissivity, r=radial distance from the well,
't=time, S=aquifer elastic storage coefficient, sigma=aquitard porosity, K"=aquitard
'permeability, B"=aquitard saturated thickness and L=length that may be chosen arbitrarily
'since it cancels out in the Boulton solution. For example, choose L=one metre if units
'of metre are used in the problem. All input and output variables are dimensionless,
'and the prime superscript has been omitted for notational convenience.
Function
W_3(r, t, epsilon, K)
Pi = 3.141592654
If
t <= 0  # Then
W_3 = 0  #
Else
If
t * epsilon * K >= 10
Then
W_3 = Exp1((1 + epsilon) * r ^ 2 / (4 * t * epsilon)) / (4 * Pi)
Else
n = 8
W_3 = 0  #
For
i = 1
To
n
W_3 = W_3 + Stehcoef(i, n) * f_3(r, i * Log(2) / t, epsilon, K)
Next
i
W_3 = W_3 * Log(2) / t
End
If
End
If
End
Function

'This computes the Laplace transform used in the calculation of drawdown for flow to
'a well in a delayed-yield aquifer.
Function
f_3(r, p, epsilon, K)
Pi = 3.141592654
M = Sqr(p + K - epsilon * K ^ 2 / (p + epsilon * K))
f_3 = BessK0(r * M) / (2 * Pi * p)
End
Function

'This calculates drawdowns in the top layer for delayed-yield flow to a well. All input
'and output variables are dimensionless and are defined As follows:
'   Eta_3=Abs(Eta)*T/Q    r' = r / L
t
'=t*T/SL^2  epsilon=S/sigma   K' = (K
"/B")L ^ 2 / T
'where Eta=drawdown, Q=well flow, T=transmissivity, r=radial distance from the well,
't=time,S=aquifer elastic storage coefficient, sigma=aquitard porosity, K"=aquitard
'permeability, B"=aquitard saturated thickness and L=characteristic horizontal
'distance. For a well in an infinite aquifer with no obvious value for L, choose
'L=sqr(T/(K"/B")). Then dimensionless program variables become
'  Eta_3=Abs(Eta)*T/Q   r' = r * sqr((K
"/B") / T)   t
'=t(K"/B")/S   epsilon=S/sigma   K' = 1
Function
Eta_3(r, t, epsilon, K)
Pi = 3.141592654
If
t <= 0  # Then
Eta_3 = 0  #
Else
If
t * epsilon * K >= 10
Then
Eta_3 = Exp1((1 + epsilon) * r ^ 2 / (4 * t * epsilon)) / (4 * Pi)
Else
n = 8
Eta_3 = 0  #
For
i = 1
To
n
Eta_3 = Eta_3 + Stehcoef(i, n) * g_3(r, i * Log(2) / t, epsilon, K)
Next
i
Eta_3 = Eta_3 * Log(2) / t
End
If
End
If
End
Function
'This computes the Laplace transform used in the calculation of aquitard free-surface
'drawdown for flow to a well in a delayed-yield aquifer.
Function
g_3(r, p, epsilon, K)
Pi = 3.141592654
M = Sqr(p + K - epsilon * K ^ 2 / (p + epsilon * K))
g_3 = epsilon * K * BessK0(r * M) / (2 * Pi * p * (p + epsilon * K))
End
Function

'Calculation of drawdowns created by a well in a delayed-yield aquifer next to a stream.
'All input and output variables are dimensionless and are defined as follows:
'       W_4=Abs(h)*T/Q   x=x/L  y=y/L  t=t*T/(S*L^2)  lambda=lambda*L/T
'                    K=(K' / B')*L^2/T  epsilon=S/sigma
'where Abs(h)=drawdown, T=transmissivity of the bottom aquifer, Q=well flow rate, L=shortest
'distance between the well and stream edge, S=storage coefficient of the bottom aquifer
'(either storativity or specific yield), x=coordinate measured from the stream edge toward the
'well (normal to the stream edge), y=coordinate measured along the stream edge,lambda=stream
'bed leakage coefficient and K" and B" = permeability and thickness of the top aquifer. All
'prime superscripts have been omitted in the program for notational convenience.
'NOTE: Setting K' = 0
for any value of epsilon gives the solution obtained by B.Hunt(1999)
'GROUND WATER,37(1),98-102 for either a completely confined or completely unconfined aquifer.
Function
W_4(x, y, t, lambda , K, epsilon)
If
t <= 0
Then
W_4 = 0  #
Else
If
lambda = 0  # Then
       W_4 = 0  #
       Else
       n = 20 'Increase accuracy by increasing n to any even integer. This will also
       'increase computational times.
       delta = 1 / n
       W_4 = 0  #
       y_1 = 0  #
       For i = 1 To n - 1 Step 2
       y_2 = integrand_3(i * delta, x, y, t, lambda, K, epsilon)
y_3 = integrand_3((i + 1) * delta, x, y, t, lambda , K, epsilon)
W_4 = W_4 + delta * (y_1 + 4 * y_2 + y_3) / 3
y_1 = y_3
Next
i
End
If
r = Sqr((x - 1) ^ 2 + y ^ 2)
W_4 = W_3(r, t, epsilon, K) - W_4
End
If
End
Function

'This computes the Laplace transform used in the calculation of the aquifer
'drawdown, W_4, for flow to a well in a delayed-yield aquifer.
Function
integrand_3(u, x, y, t, lambda , K, epsilon)
r = Sqr((1 + Abs(x) + 2 * Log(1 / u) /
lambda ) ^ 2 + y ^ 2)
       integrand_3 = W_3(r, t, epsilon, K)
End
Function

'This computes G(alpha,t)
Function
g_1(alpha, t, K, epsilon)
a = Abs(epsilon * K * t * (1 - alpha ^ 2))
b = Abs(K * t * alpha ^ 2)
g_1 = 1 - Exp(-((Sqr(a) - Sqr(b)) ^ 2)) * ExpBessI0(2 * Sqr(a * b)) - F_6(a, Sqr(b))
End
Function

'This calculates the total flow depletion lost from a stream when a well abstracts a flow
'Q from a delayed-yield (semi-confined) aquifer. All input and output variables are
'dimensionless with the following definitions:
' Q_4=flow depletion/Q  t' = t * T / (S * L ^ 2)
Lambda
'=Lambda*L/T  Epsilon=S/Sigma  K=(K' / B')*L^2/T
'NOTE: Setting K' = 0
for any value of epsilon gives the solution obtained by B.Hunt(1999)
'GROUND WATER,37(1),98-102 for either a completely confined or completely unconfined aquifer.
Function
Q_4(t, K, epsilon, lambda )
                          If t <= 0 Then
                          Q_4 = 0  #
                          Else
                          n = 8
                          Q_4 = 0  #
                          For i = 1 To n
                          Q_4 = Q_4 + Stehcoef(i, n) * g_4(i * Log(2) / t, K, epsilon, lambda )
                                                                                              Next i
                                                                                              Q_4 = Q_4 * Log(2) / t
                                                                                              End If
                                                                                              End Function

                                                                                              Function g_4(p, K,
                                                                                              epsilon, lambda )
                                                                                              m0 = Sqr(p * (p + K * (1 + epsilon)) / (p + epsilon * K))
                                                                                              g_4 = lambda * Exp(-m0) / (p * ( lambda + 2 * m0))
                                                                                              End Function

                                                                                              'Calculation of free surface drawdowns created by a well in a delayed-yield aquifer next
                                                                                              'to a stream. The free surface occurs in an aquitard, which is underlain by the pumped
                                                                                              'aquifer. All input and output variables are dimensionless and are defined as follows:
                                                                                              '       Eta_4=Abs(Eta)*T/Q   x=x/L  y=y/L  t=t*T/(S*L^2)  lambda=lambda*L/T
                                                                                              '                    K=(K' / B')*L^2/T  epsilon=S/sigma
                                                                                              'where Abs(Eta)=drawdown, T=transmissivity of the bottom aquifer, Q=well flow rate, L=shortest
                                                                                              'distance between the well and stream edge, S=storage coefficient of the bottom aquifer
                                                                                              '(either storativity or specific yield), x=coordinate measured from the stream edge toward the
                                                                                              'well (normal to the stream edge), y=coordinate measured along the stream edge,lambda=stream
                                                                                              'bed leakage coefficient and K" and B" = permeability and thickness of the top aquifer. All
                                                                                              'prime superscripts have been omitted in the program for notational convenience.
                                                                                              Function Eta_4(x, y, t,
                                                                                              lambda, K, epsilon)
If
t <= 0
Then
Eta_4 = 0  #
Else
If
lambda = 0  # Then
       Eta_4 = 0  #
       Else
       n = 20 'Increase accuracy by increasing n to any even integer. This will also
       'increase computational times.
       delta = 1 / n
       Eta_4 = 0  #
       y_1 = 0  #
       For i = 1 To n - 1 Step 2
       y_2 = f_4(i * delta, x, y, t, lambda, K, epsilon)
y_3 = f_4((i + 1) * delta, x, y, t, lambda , K, epsilon)
Eta_4 = Eta_4 + delta * (y_1 + 4 * y_2 + y_3) / 3
y_1 = y_3
Next
i
End
If
r = Sqr((x - 1) ^ 2 + y ^ 2)
Eta_4 = Eta_3(r, t, epsilon, K) - Eta_4
End
If
End
Function

'This computes the Laplace transform used in the calculation of aquitard free-surface
'drawdown Eta_4, for the stream-depletion problem.
Function
f_4(u, x, y, t, lambda , K, epsilon)
If
u <= 0
Then
f_4 = 0  #
Else
r = Sqr((1 + Abs(x) + 2 * Log(1 / u) /
lambda ) ^ 2 + y ^ 2)
       f_4 = Eta_3(r, t, epsilon, K)
End
If
End
Function

'This computes drawdown when flow is depleted from a spring. All input and output
'variables are dimensionless and are defined as follows:
'   W_5=sT/Q  x' = x / L
y
'=y/L  t' = tT / SL ^ 2
epsilon = S / sigma
K
'=(K' / B')L^2/T
'                    x_0' = x_0 / L
alpha = A(K
"/B") / T
'where s=drawdown, Q=well abstraction rate, T=pumped aquifer transmissivity,  t=time,
'S=pumped aquifer storage coefficient, L=distance between the well and spring, sigma=
'aquitard specific yield (effective porosity near the free surface), K' = aquitard
'permeability, B' = aquitard
thickness, x_0 = radius
of
the
spring
outflow
area, A = spring
'outflow area in the horizontal (x,y) plane, K"=aquitard permeability beneath the spring
'and B"=aquitard thickness beneath the spring. The spring is at the coordinate origin
'(x', y
')=(0,0), the well has the coordinates (x', y
')=(1,0) and the drawdown is calculated
'at the point (x', y
'). This uses the Stehfest algorithm to invert a Laplace transform.
'Note that alpha used in this program is the reciprocal of alpha defined in the paper
'Hunt, B. (2004). "Spring depletion solution." ASCE Journal of Hydrologic Engineering,
'9(2), 144-149.
Function
W_5(x, y, t, epsilon, K, x_0, alpha)
If
t <= 0
Then
W_5 = 0  #
Else
n = 6
W_5 = 0  #
For
i = 1
To
n
W_5 = W_5 + Stehcoef(i, n) * _
W_5_tansform(x, y, i * Log(2) / t, epsilon, K, x_0, alpha)
Next
i
W_5 = W_5 * Log(2) / t
End
If
End
Function

'This computes a numerical value for the Laplace transform of the drawdown in the
'Spring-depletion problem.
Function
W_5_tansform(x, y, p, epsilon, K, x_0, alpha)
Pi = 3.141592654
r_1 = Sqr((x - 1) ^ 2 + y ^ 2)
r_2 = Sqr(x ^ 2 + y ^ 2)
M = Sqr(p * (p + K * (1 + epsilon)) / (p + epsilon * K))
W_5_tansform = BessK0(r_1 * M) - alpha * BessK0(M) * BessK0(r_2 * M) / _
(2 * Pi + alpha * BessK0(x_0 * M))
W_5_tansform = W_5_tansform / (2 * Pi * p)
End
Function

'This computes the free surface drawdown in the aquitard when flow is depleted from a
'spring. All input and output variables are dimensionless and are defined as follows:
'   Eta_5=Eta*T/Q  x' = x / L
y
'=y/L  t' = tT / SL ^ 2
epsilon = S / sigma
K
'=(K' / B')L^2/T
'                    x_0' = x_0 / L
alpha = A(K
"/B") / T
'where Eta=drawdown, Q=well abstraction rate, T=pumped aquifer transmissivity,  t=time,
'S=pumped aquifer storage coefficient, L=distance between the well and spring, sigma=
'aquitard specific yield (effective porosity near the free surface), K' = aquitard
'permeability, B' = aquitard
thickness, x_0 = radius
of
the
spring
outflow
area, A = spring
'outflow area in the horizontal (x,y) plane, K"=aquitard permeability beneath the spring
'and B"=aquitard thickness beneath the spring. The spring is at the coordinate origin
'(x', y
')=(0,0), the well has the coordinates (x', y
')=(1,0) and the drawdown is calculated
'at the point (x', y
'). This uses the Stehfest algorithm to invert a Laplace transform.
'Note that alpha used in this program is the reciprocal of alpha defined in the paper
'Hunt, B. (2004). "Spring depletion solution." ASCE Journal of Hydrologic Engineering,
'9(2), 144-149.
Function
Eta_5(x, y, t, epsilon, K, x_0, alpha)
If
t <= 0
Then
Eta_5 = 0  #
Else
n = 6
Eta_5 = 0  #
For
i = 1
To
n
Eta_5 = Eta_5 + Stehcoef(i, n) * _
Eta_5_transform(x, y, i * Log(2) / t, epsilon, K, x_0, alpha)
Next
i
Eta_5 = Eta_5 * Log(2) / t
End
If
End
Function

'This computes a Laplace transform for Eta_5.
Function
Eta_5_transform(x, y, p, epsilon, K, x_0, alpha)
Eta_5_transform = epsilon * K * W_5_tansform(x, y, p, epsilon, K, x_0, alpha) / _
(p + epsilon * K)
End
Function

'This computes flow depleted from a spring. All input and output variables are
'dimensionless and are defined as follows:
'   Q_5=Qd(t)/Q  t' = tT / SL ^ 2
epsilon = S / sigma
K
'=(K' / B')L^2/T  x_0' = x_0 / L
'                    alpha=A(K"/B")/T
'where Qd=spring depletion, Q=well abstraction rate, t=time, T=pumped aquifer
'transmissivity, S=pumped aquifer storage coefficient, L=distance between the
'well and spring, sigma=aquitard porosity, K' = aquitard
permeability, B'=
'aquitard thickness, K"=aquitard permeability beneath the spring, B"=aquitard
'thickness beneath the spring and x_0=spring radius. This uses the Stehfest algorithm
'to invert a Laplace transform. Note that alpha used in this program is the reciprocal
'of alpha defined in the paper Hunt, B. (2004). "Spring depletion solution." ASCE Journal
'of Hydrologic Engineering, 9(2), 144-149.
Function
Q_5(t, epsilon, K, x_0, alpha)
If
t <= 0
Then
Q_5 = 0  #
Else
n = 6
Q_5 = 0  #
For
i = 1
To
n
Q_5 = Q_5 + Stehcoef(i, n) * _
Q_5_transform(i * Log(2) / t, epsilon, K, x_0, alpha)
Next
i
Q_5 = Q_5 * Log(2) / t
End
If
End
Function

'This computes a Laplace transform for Q_5.
Function
Q_5_transform(p, epsilon, K, x_0, alpha)
Pi = 3.141592654
M = Sqr(p * (p + K * (1 + epsilon)) / (p + epsilon * K))
Q_5_transform = alpha * BessK0(M) / (p * (2 * Pi + alpha * BessK0(M * x_0)))
End
Function

'This calculates the piezometric-head rise in a delayed-yield aquifer as a result of
'aquitard recharge over a rectangular area. All input and output variables are
'dimensionless and are defined as follows:
' W_6=h*T/RL^2  x' = x / L
y
'=y/L  t' = t * T / SL ^ 2
epsilon = S / sigma
K
'=(K"/B")L^2/T
'                           a' = a / L
b'=b/L
'where h=head rise, R=vertical recharge flux velocity (specific discharge),
'T=transmissivity, x and y=horizontal Cartesian coordinates, t=time, S=aquifer elastic
'storage coefficient, sigma=aquitard porosity, K"=aquitard permeability, B"=aquitard
'saturated thickness, and a and b=half-width of the rectangular recharge area in the x
'and y directions, respectively. The length L may be chosen arbitrarily unless K' = 1,
'in which case L=sqrt(T/(K"/B")). The solution for recharge to an unconfined aquifer
'with no overlying aquitard is obtained from this solution by letting epsilon approach
'infinity, in which case the solution no longer depends upon K'.Consequently, calculate
'this limiting case by letting epsilon = 10, 100, ... until the solution stops changing.
'Then S in the definition of t'
must
be
interpreted as the
effective
porosity
of
the
'unconfined aquifer.
Function
W_6(x, y, t, epsilon, K, a, b)
If
t <= 0  # Then
W_6 = 0  #
Else
delta = 0.005
y_1 = 0  #
W_6 = 0  #
For
i = 1
To
200
tau = i * delta
alpha = epsilon * K * t * Abs(1 - tau ^ 2)
beta1 = tau * Sqr(K * t)
y_2 = tau * F_6(alpha, beta1) * (Erf((a - x) / (2 * tau * Sqr(t))) + _
Erf((a + x) / (2 * tau * Sqr(t)))) *(Erf((b - y) / (2 * tau * Sqr(t)))
_
+ Erf((b + y) / (2 * tau * Sqr(t))))
W_6 = W_6 + delta * (y_1 + y_2) / 2
y_1 = y_2
Next
i
W_6 = W_6 * t / 2
End
If
End
Function

'This calculates the function F(alpha,beta) needed for the computation of piezometric head
'rises in a delayed-yield aquifer as a result of aquitard recharge.
Function
F_6(alpha, beta1)
Pi = 3.141592654
a = Sqr(alpha) - beta1
b = 2 * beta1 * Sqr(alpha)
If
16 * beta1 ^ 2 >= 1000
And
8 * b >= 1000
Then
F_6 = 1 - Erfc(a) / 2 - Exp(-a * a) / (4 * beta1 * Sqr(Pi)) - _
(Erfc(a) + a * Exp(-a * a) / Sqr(Pi)) / (16 * beta1 ^ 2)
ElseIf
beta1 = Sqr(alpha)
Then
F_6 = (1 - ExpBessI0(2 * alpha)) / 2
ElseIf
beta1 > Sqr(alpha)
Then
c = Sqr(alpha) / beta1
If
c = 0  # Then
F_6 = 0  #
Else
term1 = ExpBessI0(b)
term2 = ExpBessI1(b)
F_6 = term2 * c
n = 1
Do
n = n + 1
termn = term1 - 2 * (n - 1) * term2 / b
z = termn * c ^ n
F_6 = F_6 + z
term1 = term2
term2 = termn
Loop
Until
Abs(z) < 0.0001
F_6 = Exp(-a * a) * F_6
End
If
Else
c = beta1 / Sqr(alpha)
If
c = 0
Then
F_6 = 1 - Exp(-a * a)
Else
term1 = ExpBessI0(b)
term2 = ExpBessI1(b)
F_6 = term1 + term2 * c
n = 1
Do
n = n + 1
termn = term1 - 2 * (n - 1) * term2 / b
z = termn * c ^ n
F_6 = F_6 + z
term1 = term2
term2 = termn
Loop
Until
Abs(z) < 0.0001
F_6 = 1 - Exp(-a * a) * F_6
End
If
End
If
End
Function

'This calculates the free surface rise in an overlying aquitard for a delayed-yield
'aquifer as a result of recharge over a rectangular area. All input and output variables
'are dimensionless and are defined as follows:
' Eta_6=Eta*T/RL^2  x' = x / L
y
'=y/L  t' = t * T / SL ^ 2
epsilon = S / sigma
K
'=(K"/B")L^2/T
'                           a' = a / L
b'=b/L
'where Eta=free surface rise, R=vertical recharge flux velocity (specific discharge),
'T=transmissivity, x and y=horizontal Cartesian coordinates, t=time, S=aquifer elastic
'storage coefficient, sigma=aquitard porosity, K"=aquitard permeability, B"=aquitard
'saturated thickness, and a and b=half-width of the rectangular recharge area in the x
'and y directions, respectively. The length L may be chosen arbitrarily unless K' = 1,
'in which case L=sqrt(T/(K"/B")).
Function
Eta_6(x, y, t, epsilon, K, a, b)
If
t <= 0  # Then
Eta_6 = 0  #
Else
delta = 0.005
y_1 = 0  #
Eta_6 = 0  #
For
i = 1
To
200
tau = i * delta
alpha = epsilon * K * t * Abs(1 - tau ^ 2)
beta1 = tau * Sqr(K * t)
y_2 = tau * G_6(alpha, beta1) * (Erf((a - x) / (2 * tau * Sqr(t))) + _
Erf((a + x) / (2 * tau * Sqr(t)))) *(Erf((b - y) / (2 * tau * Sqr(t)))
_
+ Erf((b + y) / (2 * tau * Sqr(t))))
Eta_6 = Eta_6 + delta * (y_1 + y_2) / 2
y_1 = y_2
Next
i
Eta_6 = Eta_6 * t / 2
If
Abs(x) <= a
Then
If
Abs(y) <= b
Then
Eta_6 = (1 - Exp(-epsilon * K * t)) / K + Eta_6
End
If
End
If
End
If
End
Function

'This calculates the function G(alpha,beta) needed for the computation of free surface rises
'in an aquitard as a result of recharge. The aquitard overlies a delayed-yield aquifer.
Function
G_6(alpha, beta1)
a = (Sqr(alpha) - beta1) ^ 2
b = 2 * beta1 * Sqr(alpha)
If
beta1 = 0  # Then
G_6 = 1 - (alpha + 1) * Exp(-alpha)
Else
c = Sqr(alpha) / beta1
G_6 = -Exp(-a) * ExpBessI1(b) * c + F_6(alpha, beta1)
End
If
End
Function

'This calculates the piezometric-head rise in a delayed-yield aquifer as a result of
'aquitard recharge over a circular  area. All input and output variables are
'dimensionless and are defined as follows:
' W_7=h*T/RL^2  r' = r / L
t
'=t*T/SL^2  epsilon=S/sigma  K' = (K
"/B")L ^ 2 / T
a
'=a/L
'where h=head rise, R=vertical recharge flux velocity (specific discharge),
'T=transmissivity, r=radial coordinate, t=time, S=aquifer elastic storage coefficient,
'sigma=aquitard porosity, K"=aquitard permeability, B"=aquitard saturated thickness, and
'a=radius of the circular recharge area. The length L may be chosen arbitrarily unless K' = 1,
'in which case L=sqrt(T/(K"/B")). The solution for recharge to an unconfined aquifer
'with no overlying aquitard is obtained from this solution by letting epsilon approach
'infinity, in which case the solution no longer depends upon K'.Consequently, calculate
'this limiting case by letting epsilon = 10, 100, ... until the solution stops changing.
'Then S in the definition of t'
must
be
interpreted as the
effective
porosity
of
the
'unconfined aquifer.
Function
W_7(r, t, epsilon, K, a)
If
t <= 0  # Then
W_7 = 0  #
Else
delta = 0.005
y_1 = 0  #
W_7 = 0  #
For
i = 1
To
200
tau = i * delta
alpha = epsilon * K * t * Abs(1 - tau ^ 2)
beta1 = tau * Sqr(K * t)
Alpha2 = a ^ 2 / (4 * t * tau ^ 2)
beta2 = r / (2 * tau * Sqr(t))
y_2 = tau * F_6(Alpha2, beta2) * F_6(alpha, beta1)
W_7 = W_7 + delta * (y_1 + y_2) / 2
y_1 = y_2
Next
i
W_7 = 2 * t * W_7
End
If
End
Function

'This calculates the free surface rise in an overlying aquitard for a delayed-yield
'aquifer as a result of aquitard recharge over a circular  area. All input and output
'variables are dimensionless and are defined as follows:
' Eta_7=Eta*T/RL^2  r' = r / L
t
'=t*T/SL^2  epsilon=S/sigma  K' = (K
"/B")L ^ 2 / T
a
'=a/L
'where Eta=free surface rise, R=vertical recharge flux velocity (specific discharge),
'T=transmissivity, r=radial coordinate, t=time, S=aquifer elastic storage coefficient,
'sigma=aquitard porosity, K"=aquitard permeability, B"=aquitard saturated thickness, and
'a=radius of the circular recharge area. The length L may be chosen arbitrarily unless K' = 1,
'in which case L=sqrt(T/(K"/B")).
Function
Eta_7(r, t, epsilon, K, a)
If
t <= 0  # Then
Eta_7 = 0  #
Else
delta = 0.005
y_1 = 0  #
Eta_7 = 0  #
For
i = 1
To
200
tau = i * delta
alpha = epsilon * K * t * Abs(1 - tau ^ 2)
beta1 = tau * Sqr(K * t)
Alpha2 = a ^ 2 / (4 * t * tau ^ 2)
beta2 = r / (2 * tau * Sqr(t))
y_2 = tau * F_6(Alpha2, beta2) * G_6(alpha, beta1)
Eta_7 = Eta_7 + delta * (y_1 + y_2) / 2
y_1 = y_2
Next
i
Eta_7 = 2 * t * Eta_7
If
r <= a
Then
Eta_7 = (1 - Exp(-epsilon * K * t)) / K + Eta_7
End
If
End
If
End
Function

'This calculates drawdowns for flow to a vertical or non-vertical well in a leaky aquifer.
'The bottom boundary is impermeable, the top boundary is an aquitard and the well is
'approximated with a straight-line segment that has any orientation, length and location
'between the top and bottom aquifer boundaries. All input and output variables are
'dimensionless and are defined as follows:
'                              W_8=s*B*KH/Q
'  x' = x / B
y
'=y/B  z' = z / B
t
'=t*KH/Ss*B^2  K=KV/KH  lambda=(K' / B')/(KV/B)  x_0' = x_0 / B
'           y_0' = y_0 / B
z_0
'=z_0/B  x_M' = x_M / B
y_M
'=y_M/B  z_M' = z_M / B
'
'where the straight-line segment from (x_0,y_0,z_0) to (x_M,y_M,z_M) has been divided
'into M segments of equal length, B=aquifer thickness, Ss=specific storage, KH=horizontal
'permeability, KV=vertical permeability, B' = aquitard
thickness, K
'=vertical aquitard
'permeability and nterms=number of terms retained in the infinite series. A point sink has
'been placed at the midpoint of each segment, the drawdown is calculated at (x,y,z), and the
'prime superscripts are omitted in the program for notational convenience. A value of nterms=10
'is probably sufficient for most calculations. However, if unexpected oscillations appear in
'the solution, try increasing nterms and see if they disappear.
Function
W_8(x, y, z, t, K, lambda , x_0, y_0, z_0, x_M, y_M, z_M, M, nterms)
If
t <= 0  # Then
W_8 = 0  #
Else
Pi = 3.141592654
W_8 = 0  #
For
j = 1
To
M
xj = x_0 + (x_M - x_0) * (j - 1 / 2) / M
yj = y_0 + (y_M - y_0) * (j - 1 / 2) / M
zj = z_0 + (z_M - z_0) * (j - 1 / 2) / M
r = Sqr((x - xj) ^ 2 + (y - yj) ^ 2)
For
n = 1
To
nterms
a = rt(n, lambda )
                 If a = 0  # Then
                 d = 1  #
                 Else
                 d = (1 + Sin(2 * a) / (2 * a)) / 2
                 End If
                 term = w(r ^ 2 / (4 * t), r *a * Sqr(K)) * Cos(a * zj) * Cos(a * z) / (4 * Pi * d)
W_8 = W_8 + term
Next
n
Next
j
W_8 = W_8 / M
End
If
End
Function

'This uses Newton'
s
method
to
calculate
the
n
th
root
of
the
equation
'                    x*tan(x)=lambda
'where a positive value is specified for lambda, n>=1 and 0<=x<infinity.
'This routine is used in the program W_8.
Function
rt(n, lambda )
             Pi = 3.141592654
             If lambda < 0.001 Then
             If n = 1 Then
             rt = Sqr( lambda )
             Exit Function
             Else
             rt = (n - 1) * Pi + lambda / ((n - 1) * Pi)
             Exit Function
             End If
             End If
             If lambda <= 4  # Then
             x = (n - 1) * Pi + Pi / 4
             Do
             f = x * Sin(x) - lambda * Cos(x)
             fp = (1 + lambda ) * Sin(x) + x * Cos(x)
             Change = f / fp
             x = x - Change
             Loop Until Abs(Change) < 0.000001
             rt = x
             Else
             epsilon = 1 / lambda
             x = (n - 1) * Pi + Pi / 4
             Do
             f = epsilon * x * Sin(x) - Cos(x)
             fp = (1 + epsilon) * Sin(x) + epsilon * x * Cos(x)
             Change = f / fp
             x = x - Change
             Loop Until Abs(Change) < 0.000001
             rt = x
             End If
             End Function
             'This calulates the fluid pore velocity in the x direction for flow to a non-vertical well in a
             'leaky aquifer. The bottom boundary is impermeable, the top boundary is an aquitard and the
             'well is approximated with a straight-line segment that has any orientation, length and location
             'between the top and bottom aquifer boundaries. All input and output variables are dimensionless
             'and are defined as follows:
             '                       dx'dt'=Q' * Vx' [Q'=Q * Ss / (sigma * KH * B) and Vx'=Vx*B^2/Q]
             '  x'=x / B  y'=y/B  z'=z / B  t'=t*KH/Ss*B^2  K=KV/Kr  lambda=(K' / B')/(KV/B)  x_0'=x_0 / B
             '           y_0'=y_0 / B  z_0'=z_0/B  x_M'=x_M / B  y_M'=y_M/B  z_M'=z_M / B
             '
             'where the straight-line segment from (x_0,y_0,z_0) to (x_M,y_M,z_M) has been divided into M
             'segments of equal length, B=aquifer thickness, Ss=specific storage, KH=horizontal permeability,
             'KV=vertical permeability, Vx=specific discharge, B'=aquitard thickness, K'=vertical aquitard
'permeability and nterms=number of terms retained in the infinite series. A point sink has
'been placed at the midpoint of each segment, the drawdown is calculated at (x,y,z), and the
'prime superscripts are omitted in the program for notational convenience. A value of nterms=10
'is probably sufficient for most calculations. However, if unexpected oscillations appear in
'the solution, try increasing nterms and see if the disappear.
Function
dxdt(x, y, z, t, K, Q, lambda , x_0, y_0, z_0, x_M, y_M, z_M, M, nterms)
If
t = 0  # Then
dxdt = 0  #
Else
Pi = 3.141592654
dxdt = 0  #
For
j = 1
To
M
xj = x_0 + (x_M - x_0) * (j - 1 / 2) / M
yj = y_0 + (y_M - y_0) * (j - 1 / 2) / M
zj = z_0 + (z_M - z_0) * (j - 1 / 2) / M
r = Sqr((x - xj) ^ 2 + (y - yj) ^ 2)
For
n = 1
To
nterms
a = rt(n, lambda )
                 If a = 0  # Then
                 d = 1  #
                 Else
                 d = (1 + Sin(2 * a) / (2 * a)) / 2
                 End If
                 term = Cos(a * zj) * Cos(a * z) / (4 * Pi * d)
                 term = term * (Wx(r ^ 2 / (4 * t), r *a * Sqr(K)) * (x - xj) / (2 * t) + _
Wx(r ^ 2 / (4 * t), r * a * Sqr(K)) * (x - xj) * a * Sqr(K) / r)
dxdt = dxdt + term
Next
n
Next
j
dxdt = Q * dxdt / M
End
If
End
Function

'This calulates the fluid pore velocity in the y direction for flow to a non-vertical well in a
'leaky aquifer. The bottom boundary is impermeable, the top boundary is an aquitard and the
'well is approximated with a straight-line segment that has any orientation, length and location
'between the top and bottom aquifer boundaries. All input and output variables are dimensionless
'and are defined as follows:
'                       dy'
dt
'=Q' * Vy
' [Q' = Q * Ss / (sigma * KH * B) and Vy
'=Vy*B^2/Q]
'  x' = x / B
y
'=y/B  z' = z / B
t
'=t*KH/Ss*B^2  K=KV/Kr  lambda=(K' / B')/(KV/B)  x_0' = x_0 / B
'           y_0' = y_0 / B
z_0
'=z_0/B  x_M' = x_M / B
y_M
'=y_M/B  z_M' = z_M / B
'
'where the straight-line segment from (x_0,y_0,z_0) to (x_M,y_M,z_M) has been divided into M
'segments of equal length, B=aquifer thickness, Ss=specific storage, KH=horizontal permeability,
'KV=vertical permeability, Vy=specific discharge, B' = aquitard
thickness, K
'=vertical aquitard
'permeability and nterms=number of terms retained in the infinite series. A point sink has
'been placed at the midpoint of each segment, the drawdown is calculated at (x,y,z), and the
'prime superscripts are omitted in the program for notational convenience. A value of nterms=10
'is probably sufficient for most calculations. However, if unexpected oscillations appear in
'the solution, try increasing nterms and see if the disappear.
Function
dydt(x, y, z, t, K, Q, lambda , x_0, y_0, z_0, x_M, y_M, z_M, M, nterms)
If
t = 0  # Then
dydt = 0  #
Else
Pi = 3.141592654
dydt = 0  #
For
j = 1
To
M
xj = x_0 + (x_M - x_0) * (j - 1 / 2) / M
yj = y_0 + (y_M - y_0) * (j - 1 / 2) / M
zj = z_0 + (z_M - z_0) * (j - 1 / 2) / M
r = Sqr((x - xj) ^ 2 + (y - yj) ^ 2)
For
n = 1
To
nterms
a = rt(n, lambda )
                 If a = 0  # Then
                 d = 1  #
                 Else
                 d = (1 + Sin(2 * a) / (2 * a)) / 2
                 End If
                 term = Cos(a * zj) * Cos(a * z) / (4 * Pi * d)
                 term = term * (Wx(r ^ 2 / (4 * t), r *a * Sqr(K)) * (y - yj) / (2 * t) + _
Wx(r ^ 2 / (4 * t), r * a * Sqr(K)) * (y - yj) * a * Sqr(K) / r)
dydt = dydt + term
Next
n
Next
j
dydt = Q * dydt / M
End
If
End
Function

'This calulates the fluid pore velocity in the z direction for flow to a non-vertical well in a
'leaky aquifer. The bottom boundary is impermeable, the top boundary is an aquitard and the
'well is approximated with a straight-line segment that has any orientation, length and location
'between the top and bottom aquifer boundaries. All input and output variables are dimensionless
'and are defined as follows:
'                       dz'
dt
'=Q' * Vz
' [Q' = Q * Ss / (sigma * KH * B) and Vz
'=Vz*B^2/Q]
'  x' = x / B
y
'=y/B  z' = z / B
t
'=t*KH/Ss*B^2  K=KV/Kr  lambda=(K' / B')/(KV/B)  x_0' = x_0 / B
'           y_0' = y_0 / B
z_0
'=z_0/B  x_M' = x_M / B
y_M
'=y_M/B  z_M' = z_M / B
'
'where the straight-line segment from (x_0,y_0,z_0) to (x_M,y_M,z_M) has been divided into M
'segments of equal length, B=aquifer thickness, Ss=specific storage, KH=horizontal permeability,
'KV=vertical permeability, Vz=specific discharge, B' = aquitard
thickness, K
'=vertical aquitard
'permeability and nterms=number of terms retained in the infinite series. A point sink has
'been placed at the midpoint of each segment, the drawdown is calculated at (x,y,z), and the
'prime superscripts are omitted in the program for notational convenience. A value of nterms=10
'is probably sufficient for most calculations. However, if unexpected oscillations appear in
'the solution, try increasing nterms and see if the disappear.
Function
dzdt(x, y, z, t, K, Q, lambda , x_0, y_0, z_0, x_M, y_M, z_M, M, nterms)
If
t = 0  # Then
dzdt = 0  #
Else
Pi = 3.141592654
dzdt = 0  #
For
j = 1
To
M
xj = x_0 + (x_M - x_0) * (j - 1 / 2) / M
yj = y_0 + (y_M - y_0) * (j - 1 / 2) / M
zj = z_0 + (z_M - z_0) * (j - 1 / 2) / M
r = Sqr((x - xj) ^ 2 + (y - yj) ^ 2)
For
n = 1
To
nterms
a = rt(n, lambda )
                 If a = 0  # Then
                 d = 1  #
                 Else
                 d = (1 + Sin(2 * a) / (2 * a)) / 2
                 End If
                 term = -a * w(r ^ 2 / (4 * t), r *a * Sqr(K)) * Cos(a * zj) * Sin(a * z)
_
/ (4 * Pi * d)
dzdt = dzdt + term
Next
n
Next
j
dzdt = K * Q * dzdt / M
End
If
End
Function

'This calculates drawdowns for flow to either a vertical or non-vertical well in a
'Neuman-type unconfined aquifer. The bottom boundary is impermeable, the top boundary is
'a free surface and the well is approximated with a straight-line segment that has any
'orientation, length and location between the top and bottom aquifer boundaries. All input
'and output variables are dimensionless and are defined as follows:
'
'                              W_9=s*B*KH/Q
'  x' = x / B
y
'=y/B  z' = z / B
t
'=t*KH/Ss*B^2  K' = KV / KH
epsilon = (KV / KH) * (Ss * B / sigma)
x_0
'=x_0/B
'           y_0' = y_0 / B
z_0
'=z_0/B  x_M' = x_M / B
y_M
'=y_M/B  z_M' = z_M / B
'
'where the straight-line segment from (x_0,y_0,z_0) to (x_M,y_M,z_M) has been divided
'into M segments of equal length, B=aquifer thickness, Ss=specific storage, KH=horizontal
'permeability, KV=vertical permeability, B' = aquitard
thickness, K
'=vertical aquitard
'permeability and sigma=effective porosity at the free surface. A point sink has
'been placed at the midpoint of each segment, the drawdown is calculated at (x,y,z), and
'prime superscripts are omitted in the program for notational convenience.
Function
W_9(x, y, z, t, K, epsilon, x_0, y_0, z_0, x_M, y_M, z_M, M)
If
t <= 0  # Then
W_9 = 0  #
Else
Pi = 3.141592654
W_9 = 0  #
For
j = 1
To
M
xj = x_0 + (x_M - x_0) * (j - 1 / 2) / M
yj = y_0 + (y_M - y_0) * (j - 1 / 2) / M
zj = z_0 + (z_M - z_0) * (j - 1 / 2) / M
r = Sqr((x - xj) ^ 2 + (y - yj) ^ 2)
W_9 = W_9 + sink_9(r, z, t, zj, K, epsilon)
Next
j
W_9 = W_9 / M
End
If
End
Function

'This uses a numerical inversion of the Laplace transform to obtain the solution
'for flow to a point sink in a Neuman-type unconfined aquifer. A list of dimensionless
'variables follows:
'
' sink_9=s*KH*B/Q  r' = r / B
z
'=z/B  t' = t * KH / (Ss * B ^ 2)
L
'=L/B  K' = KV / KH
epsilon = (KV / KH) * (Ss * B / sigma)
'
'where s=drawdown, Q=well flow, B=aquifer thickness, KH=horizontal permeability, KV=vertical
'permeability, r=radial distance from the well, z=vertical distance above the aquifer bottom
'boundary, t=time, Ss=specific storage, L=z coordinate of the sink and sigma=effective porosity
'at the free surface. All input and output variables are dimensionless, and primes are omitted
'for notational convenience.
Function
sink_9(r, z, t, L, K, epsilon)
If
t <= 0  # Then
sink_9 = 0  #
Else
sink_9 = Neuman_inv(r, z, t, L, K, epsilon)
End
If
End
Function

'This uses the Stehfest algorithm to invert a Laplace transform. See Stehfest,H.
'(1970) "Numerical inversion of Laplace transforms," Comm. ACM, Vol.13, No.1,
'pages 47-49 and 624. The approximation has the nature of an asymptotic series,
'in which n is the number of terms used. Thus, n shuld not exceed about 20, and
'n = 10 is probably sufficient for most applications. Note that n must be an even
'integer and that t = time. Beware: this works well for some transforms but very
'poorly for others. Furthermore, a value for n that is slightly too large can
'cause a dramatic decrease in accuracy.
Function
Neuman_inv(r, z, t, L, K, epsilon)
ninv = 8
ninv = Application.Even(ninv)
Neuman_inv = 0  #
For
i = 1
To
ninv
Neuman_inv = Neuman_inv + Stehcoef(i, ninv) * _
Neuman(r, z, i * Log(2) / t, L, K, epsilon)
Next
i
Neuman_inv = Neuman_inv * Log(2) / t
End
Function

'This calculates a numerical value for the Laplace transform
'for flow to a point sink in a Neuman-type unconfined aquifer.
Function
Neuman(r, z, p, L, K, epsilon)
Pi = 3.141592654
Neuman = 0  #
n = 0
Do
n = n + 1
lambda = rt(n, p / epsilon)
er = BessK0(r * Sqr(p + K *
lambda ^ 2)) / (1 + Sin(2 * lambda ) / (2 * lambda ))
       Neuman = Neuman + Cos( lambda * L) * Cos( lambda * z) * er
       Loop Until er <= 0.000001
       Neuman = Neuman / (Pi * p)
       End Function

       'This calculates drawdowns for flow to either a vertical or non-vertical well in a
       'Boulton-type semi-confined aquifer. A pumped aquifer is overlain by an aquitard containing
       'a standing water table, and the well is approximated with a straight-line segment that has any
       'orientation, length and location between the top and bottom aquifer boundaries. The free surface
       'in the aquitard will drawd down with time. All input and output variables are dimensionless
       'and are defined as follows:
       '
       '                              W_10=s*B*KH/Q
       '  x'=x / B  y'=y/B  z'=z / B  t'=t*KH/Ss*B^2  K'=KV / KH  a=(KV / B) / (KV'/B')
       '            epsilon=((KV' / B')/(KH/B))*(Ss*B/sigma)
       ' x_0'=x_0 / B  y_0'=y_0/B  z_0'=z_0 / B  x_M'=x_M/B  y_M'=y_M / B  z_M'=z_M/B
       '
       'where the straight-line segment from (x_0,y_0,z_0) to (x_M,y_M,z_M) has been divided
       'into M segments of equal length, B=aquifer thickness, Ss=specific storage, KH=horizontal aquifer
       'permeability, KV=vertical aquifer permeability, B'=aquitard thickness, KV'=vertical aquitard
'permeability and sigma=effective porosity of the aquitard at the free surface. A point sink has
'been placed at the midpoint of each segment, the drawdown is calculated at (x,y,z), and
'prime superscripts are omitted in the program for notational convenience.
Function
W_10(x, y, z, t, K, a, epsilon, x_0, y_0, z_0, x_M, y_M, z_M, M)
If
t <= 0  # Then
W_10 = 0  #
Else
Pi = 3.141592654
W_10 = 0  #
For
j = 1
To
M
xj = x_0 + (x_M - x_0) * (j - 1 / 2) / M
yj = y_0 + (y_M - y_0) * (j - 1 / 2) / M
zj = z_0 + (z_M - z_0) * (j - 1 / 2) / M
r = Sqr((x - xj) ^ 2 + (y - yj) ^ 2)
W_10 = W_10 + sink_10(r, z, t, zj, K, a, epsilon)
Next
j
W_10 = W_10 / M
End
If
End
Function

'This uses a numerical inversion of the Laplace transform to obtain the solution
'for flow to a point sink in a Boulton-type semi-confined aquifer. A list of dimensionless
'variables follows:
'
'                              sink_10=s*B*KH/Q
'     r' = r / B
z
'=z/B  t' = t * KH / Ss * B ^ 2
L
'=L/B  K' = KV / KH
a = (KV / B) / (KV
'/B')
'            epsilon=((KV' / B')/(KH/B))*(Ss*B/sigma)
'
'where s=drawdown, Q=well flow, B=aquifer thickness, KH=horizontal permeability, KV=vertical
'permeability, r=radial distance from the well, z=vertical distance above the aquifer bottom
'boundary, t=time, Ss=specific storage, L=z coordinate of the sink and sigma=effective porosity
'at the free surface. All input and output variables are dimensionless, and primes are omitted
'for notational convenience.
Function
sink_10(r, z, t, L, K, a, epsilon)
If
t <= 0  # Then
sink_10 = 0  #
Else
sink_10 = Boulton_inv(r, z, t, L, K, a, epsilon)
End
If
End
Function

'This uses the Stehfest algorithm to invert a Laplace transform. See Stehfest,H.
'(1970) "Numerical inversion of Laplace transforms," Comm. ACM, Vol.13, No.1,
'pages 47-49 and 624. The approximation has the nature of an asymptotic series,
'in which n is the number of terms used. Thus, n shuld not exceed about 20, and
'n = 10 is probably sufficient for most applications. Note that n must be an even
'integer and that t = time. Beware: this works well for some transforms but very
'poorly for others. Furthermore, a value for n that is slightly too large can
'cause a dramatic decrease in accuracy.
Function
Boulton_inv(r, z, t, L, K, a, epsilon)
ninv = 6
ninv = Application.Even(ninv)
Boulton_inv = 0  #
For
i = 1
To
ninv
Boulton_inv = Boulton_inv + Stehcoef(i, ninv) * _
Boulton(r, z, i * Log(2) / t, L, K, a, epsilon)
Next
i
Boulton_inv = Boulton_inv * Log(2) / t
End
Function

'This calculates drawdowns of the free surface within an aquitard for flow to either a vertical
'or non-vertical well screened in an underlying aquifer. This type of aquifer is referred to as
'a Boulton-type semi-confined aquifer. A pumped aquifer is overlain by an aquitard containing
'a standing water table, and the well is approximated with a straight-line segment that has any
'orientation, length and location between the top and bottom aquifer boundaries. The free surface
'in the aquitard will drawd down with time. All input and output variables are dimensionless
'and are defined as follows:
'
'                              Eta_10=Eta*B*KH/Q
'  x' = x / B
y
'=y/B  z' = z / B
t
'=t*KH/Ss*B^2  K' = KV / KH
a = (KV / B) / (KV
'/B')
'            epsilon=((KV' / B')/(KH/B))*(Ss*B/sigma)
' x_0' = x_0 / B
y_0
'=y_0/B  z_0' = z_0 / B
x_M
'=x_M/B  y_M' = y_M / B
z_M
'=z_M/B
'
'where the straight-line segment from (x_0,y_0,z_0) to (x_M,y_M,z_M) has been divided
'into M segments of equal length, B=aquifer thickness, Ss=specific storage, KH=horizontal aquifer
'permeability, KV=vertical aquifer permeability, B' = aquitard
thickness, KV
'=vertical aquitard
'permeability and sigma=effective porosity of the aquitard at the free surface. A point sink has
'been placed at the midpoint of each segment, the drawdown is calculated at (x,y,z), and
'prime superscripts are omitted in the program for notational convenience.
Function
Eta_10(x, y, z, t, K, a, epsilon, x_0, y_0, z_0, x_M, y_M, z_M, M)
If
t <= 0  # Then
Eta_10 = 0  #
Else
Pi = 3.141592654
Eta_10 = 0  #
For
j = 1
To
M
xj = x_0 + (x_M - x_0) * (j - 1 / 2) / M
yj = y_0 + (y_M - y_0) * (j - 1 / 2) / M
zj = z_0 + (z_M - z_0) * (j - 1 / 2) / M
r = Sqr((x - xj) ^ 2 + (y - yj) ^ 2)
Eta_10 = Eta_10 + fs_10(r, z, t, zj, K, a, epsilon)
Next
j
Eta_10 = Eta_10 / M
End
If
End
Function

'This uses a numerical inversion of the Laplace transform to obtain the free-surface solution
'for flow to a point sink in a Boulton-type semi-confined aquifer. A list of dimensionless
'variables follows:
'
'                              fs_10=Eta*B*KH/Q
'     r' = r / B
z
'=z/B  t' = t * KH / Ss * B ^ 2
L
'=L/B  K' = KV / KH
a = (KV / B) / (KV
'/B')
'            epsilon=((KV' / B')/(KH/B))*(Ss*B/sigma)
'
'where s=drawdown, Q=well flow, B=aquifer thickness, KH=horizontal permeability, KV=vertical
'permeability, r=radial distance from the well, z=vertical distance above the aquifer bottom
'boundary, t=time, Ss=specific storage, L=z coordinate of the sink and sigma=effective porosity
'at the free surface. All input and output variables are dimensionless, and primes are omitted
'for notational convenience.
Function
fs_10(r, z, t, L, K, a, epsilon)
If
t <= 0  # Then
fs_10 = 0  #
Else
fs_10 = Boulton_Eta_inv(r, t, L, K, a, epsilon)
End
If
End
Function

'This uses the Stehfest algorithm to invert a Laplace transform. See Stehfest,H.
'(1970) "Numerical inversion of Laplace transforms," Comm. ACM, Vol.13, No.1,
'pages 47-49 and 624. The approximation has the nature of an asymptotic series,
'in which n is the number of terms used. Thus, n shuld not exceed about 20, and
'n = 10 is probably sufficient for most applications. Note that n must be an even
'integer and that t = time. Beware: this works well for some transforms but very
'poorly for others. Furthermore, a value for n that is slightly too large can
'cause a dramatic decrease in accuracy.
Function
Boulton_Eta_inv(r, t, L, K, a, epsilon)
ninv = 6
ninv = Application.Even(ninv)
Boulton_Eta_inv = 0  #
For
i = 1
To
ninv
Boulton_Eta_inv = Boulton_Eta_inv + Stehcoef(i, ninv) * _
Boulton_Eta(r, i * Log(2) / t, L, K, a, epsilon)
Next
i
Boulton_Eta_inv = Boulton_Eta_inv * Log(2) / t
End
Function

'This calculates a numerical value for the Laplace transform
'for free-surface drawdowns in flow to a point sink in a Boulton-type
'semi-confined aquifer.
Function
Boulton_Eta(r, p, L, K, a, epsilon)
Pi = 3.141592654
Boulton_Eta = 0  #
n = 0
Do
n = n + 1
d = (p / a) / (p + epsilon)
lambda = rt(n, d)
er = BessK0(r * Sqr(p + K *
lambda ^ 2)) / (1 + Sin(2 * lambda ) / (2 * lambda ))
       Boulton_Eta = Boulton_Eta + Cos( lambda * L) * er
       Loop Until er <= 0.000001
       Boulton_Eta = Boulton_Eta * epsilon / (Pi * p * (p + epsilon))
       End Function

       'This calculates a numerical value for the Laplace transform
       'for flow to a point sink in a Boulton-type semi-confined aquifer.
       Function Boulton(r, z, p, L, K, a, epsilon)
Pi = 3.141592654
d = (p / a) / (p + epsilon)
Boulton = 0  #
n = 0
Do
n = n + 1
lambda = rt(n, d)
er = BessK0(r * Sqr(p + K *
lambda ^ 2)) / (1 + Sin(2 * lambda ) / (2 * lambda ))
       Boulton = Boulton + Cos( lambda * L) * Cos( lambda * z) * er
       Loop Until er <= 0.000001
       Boulton = Boulton / (Pi * p)
       End Function

       'This calculates drawdowns in a pumped aquifer bounded above by an aquitard when the aquitard
       'lies beneath an unpumped unconfined aquifer. Dimensionless variables are defined as follows:
       '
       ' W_11=Abs(h)*T/Q  r'=r / L  t'=t*T/SL^2  epsilon=S/sigma  K'=(K"/B")L ^ 2 / T  T0'=T0/T
       '
       'where Abs(h)=drawdown, Q=well flow, T=pumped aquifer transmissivity, r=radial distance from
       'the well, L=arbitrarily chosen length, t=time, S=pumped aquifer elastic storage coefficient,
       'sigma=top unconfined aquifer porosity, K"=aquitard permeability, B"=aquitard saturated
       'thickness and T0=unconfined aquifer transmissivity. All input and output variables are
       'dimensionless, and the prime superscript has been omitted for notational convenience.
       'n = number of terms used in the Stehfest inversion algorithm. n must be an even number.
       Function W_11(r, t, epsilon, K, T0)
If
t <= 0
Then
W_11 = 0  #
Else
n = 8
W_11 = 0  #
For
i = 1
To
n
W_11 = W_11 + Stehcoef(i, n) * f_11(r, i * Log(2) / t, epsilon, K, T0)
Next
i
W_11 = W_11 * Log(2) / t
End
If
End
Function

'This calculates the Laplace transform of W_11.
Function
f_11(r, p, epsilon, K, T0)
Pi = 3.141592654
b = (p + K + (p / epsilon + K) / T0) / 2
c = (p + K * (1 + epsilon)) * p / (epsilon * T0)
lam1 = b + Sqr(b * b - c)
lam2 = b - Sqr(b * b - c)
v1 = (p + K - lam1) / K
v2 = (p + K - lam2) / K
L1 = Sqr(1 + T0 * v1 ^ 2)
L2 = Sqr(1 + T0 * v2 ^ 2)
f_11 = BessK0(r * Sqr(lam1)) / L1 ^ 2 + BessK0(r * Sqr(lam2)) / L2 ^ 2
f_11 = f_11 / (2 * Pi * p)
End
Function

'This calculates the accumulated volume of water abstracted (up to any specified time) from
'the bottom aquifer when pumping from a well screened in the bottom aquifer of a two-aquifer
'system. The dimensionless volume is defined as
'
'                            W_Volume_11=Volume*T/Sigma*Q*L^2
'
'where all dimensional variables have been defined previously in W_11.
Function
W_Volume_11(t, epsilon, K, T0)
W_Volume_11 = epsilon * t - Eta_Volume_11(t, epsilon, K, T0)
End
Function

'This calculates drawdowns in an unconfined aquifer bounded below by an aquitard when a pumped
'aquifer lies beneath the aquitard. Dimensionless variables are defined as follows:
'
' Eta_11=Abs(Eta)*T/Q  r' = r / L
t
'=t*T/SL^2  epsilon=S/sigma  K' = (K
"/B")L ^ 2 / T
T0
'=T0/T
'
'where Abs(Eta)=drawdown, Q=well flow, T=pumped aquifer transmissivity, r=radial distance from
'the well, L=arbitrarily chosen length, t=time, S=pumped aquifer elastic storage coefficient,
'sigma=top unconfined aquifer porosity, K"=aquitard permeability, B"=aquitard saturated
'thickness and T0=unconfined aquifer transmissivity. All input and output variables are
'dimensionless, and the prime superscript has been omitted for notational convenience.
'n = number of terms used in the Stehfest inversion algorithm. n must be an even number.
Function
Eta_11(r, t, epsilon, K, T0)
If
t <= 0
Then
Eta_11 = 0  #
Else
n = 8
Eta_11 = 0  #
For
i = 1
To
n
Eta_11 = Eta_11 + Stehcoef(i, n) * g_11(r, i * Log(2) / t, epsilon, K, T0)
Next
i
Eta_11 = Eta_11 * Log(2) / t
End
If
End
Function

'This calculates an integrand for the calculation of Eta_11.
Function
g_11(r, p, epsilon, K, T0)
Pi = 3.141592654
b = (p + K + (p / epsilon + K) / T0) / 2
c = (p + K * (1 + epsilon)) * p / (epsilon * T0)
lam1 = b + Sqr(b * b - c)
lam2 = b - Sqr(b * b - c)
v1 = (p + K - lam1) / K
v2 = (p + K - lam2) / K
L1 = Sqr(1 + T0 * v1 ^ 2)
L2 = Sqr(1 + T0 * v2 ^ 2)
g_11 = BessK0(r * Sqr(lam1)) * v1 / L1 ^ 2 + BessK0(r * Sqr(lam2)) * v2 / L2 ^ 2
g_11 = g_11 / (2 * Pi * p)
End
Function

'This calculates the accumulated volume of water abstracted (up to any specified time) from
'the top aquifer when pumping from a well screened in the bottom aquifer of a two-aquifer
'system. The dimensionless volume is defined as
'
'                            Eta_Volume_11=Volume*T/Sigma*Q*L^2
'
'where all dimensional variables have been defined previously in W_11.
Function
Eta_Volume_11(t, epsilon, K, T0)
If
t <= 0
Then
Eta_Volume_11 = 0  #
Else
n = 8
Eta_Volume_11 = 0  #
For
i = 1
To
n
Eta_Volume_11 = Eta_Volume_11 + Stehcoef(i, n) * _
Vol_11(i * Log(2) / t, epsilon, K, T0)
Next
i
Eta_Volume_11 = Eta_Volume_11 * Log(2) / t
End
If
End
Function

'This calculates the Laplace transform of Eta_Volume_11.
Function
Vol_11(p, epsilon, K, T0)
Pi = 3.141592654
b = (p + K + (p / epsilon + K) / T0) / 2
c = (p + K * (1 + epsilon)) * p / (epsilon * T0)
lam1 = b + Sqr(b * b - c)
lam2 = b - Sqr(b * b - c)
v1 = (p + K - lam1) / K
v2 = (p + K - lam2) / K
L1 = Sqr(1 + T0 * v1 ^ 2)
L2 = Sqr(1 + T0 * v2 ^ 2)
Vol_11 = v1 / (p * lam1 * L1 ^ 2) + v2 / (p * lam2 * L2 ^ 2)
End
Function

'This calculates drawdowns in a pumped aquifer bounded above by an aquitard when the aquitard
'lies beneath an unpumped unconfined aquifer. Flow in the aquitard is assumed compressible.
'Dimensionless variables are defined as follows:
'
' W_12=Abs(h)*T/Q  r' = r / L
t
'=t*T/SL^2  epsilon=S/sigma  K' = (K
"/B")L ^ 2 / T
T0
'=T0/T
'                                delta=B"Ss"/S
'
'where Abs(h)=drawdown, Q=well flow, T=pumped aquifer transmissivity, r=radial distance from
'the well, L=arbitrarily chosen length, t=time, S=pumped aquifer elastic storage coefficient,
'sigma=top unconfined aquifer porosity, K"=aquitard permeability, B"=aquitard saturated
'thickness, Ss"=aquitard specific storage and T0=unconfined aquifer transmissivity. All input
'and output variables are dimensionless, and the prime superscript has been omitted for
'notational convenience. n = number of terms used in the Stehfest inversion algorithm.
'n must be an even number.
Function
W_12(r, t, epsilon, K, T0, delta)
If
t <= 0
Then
W_12 = 0  #
Else
n = 8
W_12 = 0  #
For
i = 1
To
n
W_12 = W_12 + Stehcoef(i, n) * f_12(r, i * Log(2) / t, epsilon, K, T0, delta)
Next
i
W_12 = W_12 * Log(2) / t
End
If
End
Function

'This calculates the Laplace transform of W_12.
Function
f_12(r, p, epsilon, K, T0, delta)
Pi = 3.141592654
If
p * delta / K >= 100
Then
f_12 = BessK0(r * Sqr(p + Sqr(p * K * delta))) / (2 * Pi * p)
Else
a = Sqr(p * delta / K)
If
a = 0  # Then
a2 = K
a1 = a2
Else
a2 = a * K / Application.Sinh(a)
a1 = a2 * Application.Cosh(a)
End
If
b = ((p + a1) + (p / epsilon + a1) / T0) / 2
c = ((p + a1) * (p / epsilon + a1) - a2 ^ 2) / T0
lam1 = b + Sqr(b * b - c)
lam2 = b - Sqr(b * b - c)
v1 = (lam1 - (p + a1)) / a2
v2 = (lam2 - (p + a1)) / a2
L1 = Sqr(1 + T0 * v1 ^ 2)
L2 = Sqr(1 + T0 * v2 ^ 2)
f_12 = BessK0(r * Sqr(lam1)) / L1 ^ 2 + BessK0(r * Sqr(lam2)) / L2 ^ 2
f_12 = f_12 / (2 * Pi * p)
End
If
End
Function

'This calculates drawdowns in an unconfined aquifer bounded below by an aquitard when a pumped
'aquifer lies beneath the aquitard. Flow in the aquitard is assumed compressible.
'Dimensionless variables are defined as follows:
'
' Eta_12=Abs(Eta)*T/Q  r' = r / L
t
'=t*T/SL^2  epsilon=S/sigma  K' = (K
"/B")L ^ 2 / T
T0
'=T0/T
'                                  delta=B"Ss"/S
'
'where Abs(Eta)=drawdown, Q=well flow, T=pumped aquifer transmissivity, r=radial distance from
'the well, L=arbitrarily chosen length, t=time, S=pumped aquifer elastic storage coefficient,
'sigma=top unconfined aquifer porosity, K"=aquitard permeability, B"=aquitard saturated
'thickness, Ss"=aquitard specific storage and T0=unconfined aquifer transmissivity. All input
'and output variables are dimensionless, and the prime superscript has been omitted for
'notational convenience. n = number of terms used in the Stehfest inversion algorithm.
'n must be an even number.
Function
Eta_12(r, t, epsilon, K, T0, delta)
If
t <= 0
Then
Eta_12 = 0  #
Else
n = 8
Eta_12 = 0  #
For
i = 1
To
n
Eta_12 = Eta_12 + Stehcoef(i, n) * g_12(r, i * Log(2) / t, epsilon, K, T0, delta)
Next
i
Eta_12 = Eta_12 * Log(2) / t
End
If
End
Function

'This calculates an integrand for the calculation of Eta_12.
Function
g_12(r, p, epsilon, K, T0, delta)
Pi = 3.141592654
If
p * delta / K >= 100
Then
g_12 = 0  #
Else
a = Sqr(p * delta / K)
If
a = 0  # Then
a2 = K
a1 = a2
Else
a2 = a * K / Application.Sinh(a)
a1 = a2 * Application.Cosh(a)
End
If
b = ((p + a1) + (p / epsilon + a1) / T0) / 2
c = ((p + a1) * (p / epsilon + a1) - a2 ^ 2) / T0
lam1 = b + Sqr(b * b - c)
lam2 = b - Sqr(b * b - c)
v1 = (-lam1 + (p + a1)) / a2
v2 = (-lam2 + (p + a1)) / a2
L1 = Sqr(1 + T0 * v1 ^ 2)
L2 = Sqr(1 + T0 * v2 ^ 2)
g_12 = BessK0(r * Sqr(lam1)) * v1 / L1 ^ 2 + BessK0(r * Sqr(lam2)) * v2 / L2 ^ 2
g_12 = g_12 / (2 * Pi * p)
End
If
End
Function

'Calculation of drawdowns created by a well in a delayed-yield aquifer next to a stream
'when the pumped aquifer and stream both have a finite width. All input and output variables
'are dimensionless with the following definitions:
' W_13=Abs(h)*T/Q     x' = x / L
y
'=y/L        t' = t * T / (S * L ^ 2)
Epsilon = S / Sigma
'                       K=(K' / B')*L^2/T   Ks=(K"/B")*L^2/T
'where Abs(h)=drawdown, T=transmissivity of the bottom aquifer, Q=well flow rate, t=time,
'L=shortest distance between the well and stream edge, S=storage coefficient of the bottom
'aquifer (storativity), sigma=specific yield of the aquitard overlying the pumped aquifer,
'x=coordinate measured from the stream edge nearest to the well and positive toward the
'well, y=coordinate measured along the stream edge, K' / B'=permeability/thickness of the
'top aquitard not beneath the stream and K"/B"=permeability/thickness of the top aquitard
'directly beneath the stream. The portion of the pumped aquifer containing the well occupies
'the region 0<x<alpha*L for -infinity<y<infinity, and this same aquifer occupies the region
'-beta*L<x<-gama*L on the opposite side of the stream. The portion of the pumped aquifer
'beneath the stream occupies the region -gama*L<x<0. Therefore, -beta*L<-gama*L<0<alpha*L.
'All prime superscripts have been omitted in the program for notational convenience.
Function
W_13(x, y, t, K, Ks, epsilon, alpha, beta, gama)
If
t <= 0
Then
W_13 = 0  #
Else
n = 6
W_13 = 0  #
For
i = 1
To
n
W_13 = W_13 + Stehcoef(i, n) * transform_13(x, y, i * Log(2) / t, K, Ks, epsilon, _
alpha, beta, gama)
Next
i
W_13 = W_13 * Log(2) / t
End
If
End
Function

'This evaluates a Laplace transform for W_13.
Function
transform_13(x, y, p, K, Ks, epsilon, alpha, beta, gama)
Pi = 3.141592654
n = 100
umax = 20
delta = umax / n
xi = Abs(y) * delta
If
xi > 0
Then
c0 = (Sin(xi) / xi - Cos(xi)) * 4 / xi ^ 2
c1 = (Cos(xi) - Sin(xi) / xi) / xi
c2 = ((1 - 2 / xi ^ 2) * Sin(xi) + 2 * Cos(xi) / xi) / xi
Else
c0 = 4 / 3
c1 = 0  #
c2 = 1 / 3
End
If
Integral_13 = 0  #
y_1 = g_13(0, x, p, K, Ks, epsilon, alpha, beta, gama)
For
i = 1
To
n - 1
Step
2
a_13 = c2 * Cos(i * xi) - c1 * Sin(i * xi)
b_13 = c0 * Cos(i * xi)
c_13 = c1 * Sin(i * xi) + c2 * Cos(i * xi)
y_2 = g_13(i * delta, x, p, K, Ks, epsilon, alpha, beta, gama)
y_3 = g_13((i + 1) * delta, x, p, K, Ks, epsilon, alpha, beta, gama)
Integral_13 = Integral_13 + delta * (a_13 * y_1 + b_13 * y_2 + c_13 * y_3)
y_1 = y_3
Next
i
transform_13 = Integral_13
End
Function

'This evaluates an integrand for Integral_13.
Function
g_13(u, x, p, K, Ks, epsilon, alpha, beta, gama)
Pi = 3.141592654
m0 = Sqr(p * (p + K * (epsilon + 1)) / (p + epsilon * K))
M = Sqr(u ^ 2 + m0 ^ 2)
ms = Sqr(u ^ 2 + p + Ks)
denom = 2 * Pi * p * M * (esinh(M * alpha) * esinh(M * (beta - gama)) * esinh(gama * ms)
_
+ (ms / M) * ecosh(gama * ms) * esinh(M * (alpha + beta - gama))
_
+ ecosh(M * alpha) * ecosh(M * (beta - gama)) * esinh(ms * gama) * (ms / M) ^ 2)
If
x >= 1
Then
a = ecosh(M) * esinh(M * (beta - gama)) * esinh(ms * gama)
_
+ (ms / M) * ecosh(M * (1 + beta - gama)) * ecosh(ms * gama)
_
+ ecosh(M * (beta - gama)) * esinh(M) * esinh(ms * gama) * (ms / M) ^ 2
g_13 = a * ecosh(M * (alpha - x)) * Exp(-M * (x - 1)) / denom
ElseIf
x > 0
And
x < 1
Then
b = (esinh(M * (beta - gama)) * esinh(ms * gama)
_
+ (ms / M) * ecosh(M * (beta - gama)) * ecosh(ms * gama)) *ecosh(M * (alpha - 1))
c = (ecosh(ms * gama) * esinh(M * (beta - gama))
_
+ (ms / M) * ecosh(M * (beta - gama)) * esinh(ms * gama)) *(ms / M)
_
*ecosh(M * (alpha - 1))
g_13 = (b * ecosh(M * x) + c * esinh(M * x)) * Exp(-M * (1 - x)) / denom
ElseIf
x <= 0
And
x > -gama
Then
d = (esinh(M * (beta - gama)) * esinh(ms * gama)
_
+ (ms / M) * ecosh(M * (beta - gama)) * ecosh(ms * gama)) *ecosh(M * (alpha - 1))
f = (ecosh(ms * gama) * esinh(M * (beta - gama))
_
+ (ms / M) * ecosh(M * (beta - gama)) * esinh(ms * gama)) *ecosh(M * (alpha - 1))
g_13 = (d * ecosh(ms * x) + f * esinh(ms * x)) * Exp(-M + ms * x) / denom
Else
g = (ms / M) * ecosh(M * (alpha - 1))
g_13 = g * ecosh(M * (beta + x)) * Exp(-ms * gama - M * (1 - x - gama)) / denom
End
If
g_13 = 2 * g_13
End
Function

'Calculation of free surface drawdowns in an aquitard created by a well in a delayed-yield
'aquifer next to a stream when the pumped aquifer and stream both have a finite width. All
'input and output variables are dimensionless with the following definitions:
'  Eta_13=Abs(eta)*T/Q   x' = x / L
y
'=y/L     t' = t * T / (S * L ^ 2)
Epsilon = S / Sigma
'                       K=(K' / B')*L^2/T   Ks=(K"/B")*L^2/T
'where Abs(eta)=drawdown, T=transmissivity of the bottom aquifer, Q=well flow rate, t=time,
'L=shortest distance between the well and stream edge, S=storage coefficient of the bottom
'aquifer (storativity), sigma=specific yield of the aquitard overlying the pumped aquifer,
'x=coordinate measured from the stream edge nearest to the well and positive toward the
'well, y=coordinate measured along the stream edge, K' / B'=permeability/thickness of the
'top aquitard away from the stream and K"/B"=permeability/thickness of the top aquitard
'directly beneath the stream. The portion of the pumped aquifer containing the well occupies
'the region 0<x<alpha*L for -infinity<y<infinity, and this same aquifer occupies the region
'-beta*L<x<-gama*L on the opposite side of the stream. The portion of the pumped aquifer
'beneath the stream occupies the region -gama*L<x<0. Therefore, -beta*L<-gama*L<0<alpha*L.
'All prime superscripts have been omitted in the program for notational convenience.
Function
Eta_13(x, y, t, K, Ks, epsilon, alpha, beta, gama)
If
t <= 0
Then
Eta_13 = 0  #
Else
n = 6
Eta_13 = 0  #
For
i = 1
To
n
Eta_13 = Eta_13 + Stehcoef(i, n) * transform_eta_13(x, y, i * Log(2) / t, K, Ks, _
epsilon, alpha, beta, gama)
Next
i
Eta_13 = Eta_13 * Log(2) / t
End
If
End
Function

'This evaluates a Laplace transform for Eta_13.
Function
transform_eta_13(x, y, p, K, Ks, epsilon, alpha, beta, gama)
transform_eta_13 = transform_13(x, y, p, K, Ks, epsilon, alpha, beta, gama) * epsilon * K / _
(p + epsilon * K)
End
Function

'This calculates the total flow depletion lost from a stream when a well abstracts a flow
'Q from a delayed-yield (semi-confined) aquifer when the pumped aquifer has a finite width.
'All input and output variables are dimensionless with the following definitions:
'Q_13=flow depletion/Q  t' = t * T / (S * L ^ 2)
K = (K
'/B')*L ^ 2 / T
Ks = (K
"/B")*L ^ 2 / T
Epsilon = S / Sigma
'where T=transmissivity of the bottom aquifer, Q=well flow rate, t=time, L=shortest
'distance between the well and stream edge, S=storage coefficient of the bottom aquifer
'(storativity), sigma=specific yield of the aquitard overlying the pumped aquifer, K' and
'B' = permeability and thickness
of
the
aquitard
at
points
not beneath
the
stream and
'K" and B" = permeability and thickness of the aquitard beneath the stream. The aquifer
'portion containing the well occupies the region 0<x<alpha*L for -infinity<y<infinity.
'This same aquifer occupies the region -beta*L<x<-gama*L on the opposite side of the stream,
'and the stream occupies the region -gama*L<x<0 All prime superscripts have been omitted in
'the program for notational convenience.
Function
Q_13(t, K, Ks, epsilon, alpha, beta, gama)
If
t <= 0
Then
Q_13 = 0  #
Else
n = 6
Q_13 = 0  #
For
i = 1
To
n
Q_13 = Q_13 + Stehcoef(i, n) * f_13(i * Log(2) / t, K, Ks, epsilon, alpha, beta, gama)
Next
i
Q_13 = Q_13 * Log(2) / t
End
If
End
Function

'This computes a Laplace transform for Q_13, where p=transform parameter.
Function
f_13(p, K, Ks, epsilon, alpha, beta, gama)
Pi = 3.141592654
u = 0  #
m0 = Sqr(p * (p + K * (epsilon + 1)) / (p + epsilon * K))
M = Sqr(u ^ 2 + m0 ^ 2)
ms = Sqr(u ^ 2 + p + Ks)
denom = 2 * Pi * p * M * (esinh(M * alpha) * esinh(M * (beta - gama)) * esinh(gama * ms)
_
+ (ms / M) * ecosh(gama * ms) * esinh(M * (alpha + beta - gama))
_
+ ecosh(M * alpha) * ecosh(M * (beta - gama)) * esinh(ms * gama) * (ms / M) ^ 2)
d = (esinh(M * (beta - gama)) * esinh(ms * gama)
_
+ (ms / M) * ecosh(M * (beta - gama)) * ecosh(ms * gama)) *ecosh(M * (alpha - 1))
f = (ecosh(ms * gama) * esinh(M * (beta - gama))
_
+ (ms / M) * ecosh(M * (beta - gama)) * esinh(ms * gama)) *ecosh(M * (alpha - 1))
f_13 = 2 * Pi * (Ks / ms) * Exp(-M) * (d * Exp(ms * gama) * esinh(ms * gama)
_
+ f * (1 - Exp(ms * gama) * ecosh(ms * gama))) / denom
End
Function

'This computes the drawdowns in the pumped aquifer for flow to a well in a Boulton-type
'delayed-yield aquifer that occupies the region -beta*L<x<alpha*L for -infinity<y<infinity
'in the horizontal (x, y) plane. The well is located at (x, y)=(0,0). The pumped aquifer
'is overlain by an aquitard containing a free surface. All input and output variables are
'dimensionless and are defined as follows:
'  W_14=sT/Q   x' = x / L
y
'=y/L   t' = tT / SL ^ 2
K = (K
'/B')L ^ 2 / T
epsilon = S / sigma
'where s=drawdown, T=transmissivity and S=storativity of the pumped aquifer. Also, t=time,
'K' = aquitard
hydraulic
conductivity, B'=aquitard saturated thickness, L=length chosen
'arbitrarily and sigma=aquitard specific yield or effective porosity at the location of the
'free surface. This uses the method of images, and nmax controls the number of image points
'used. (Increasing nmax by one adds four new image points to the image well system.) The
'prime superscript has been omitted for notational convenience.
Function
W_14(x, y, t, K, epsilon, alpha, beta, nmax)
If
t <= 0  # Then
W_14 = 0  #
Else
r = Sqr(x ^ 2 + y ^ 2)
W_14 = W_3(r, t, epsilon, K)
For
n = 0
To
nmax
x1 = alpha - (x + beta) + (2 * n + 1) * (alpha + beta)
x2 = alpha + (x + beta) + (2 * n + 1) * (alpha + beta)
x3 = beta - (alpha - x) + (2 * n + 1) * (alpha + beta)
x4 = beta + (alpha - x) + (2 * n + 1) * (alpha + beta)
r1 = Sqr(x1 ^ 2 + y ^ 2)
r2 = Sqr(x2 ^ 2 + y ^ 2)
r3 = Sqr(x3 ^ 2 + y ^ 2)
r4 = Sqr(x4 ^ 2 + y ^ 2)
W_14 = W_14 + W_3(r1, t, epsilon, K) + W_3(r2, t, epsilon, K) + _
W_3(r3, t, epsilon, K) + W_3(r4, t, epsilon, K)
Next
n
End
If
End
Function

'This computes free surface drawdowns in the aquitard for flow to a well in a Boulton-type
'delayed-yield aquifer that occupies the region -beta*L<x<alpha*L for -infinity<y<infinity
'in the horizontal (x, y) plane. The well is located at (x, y)=(L,0). The pumped aquifer
'is overlain by an aquitard containing a free surface. All input and output variables are
'dimensionless and are defined as follows:
'  Eta_14=Eta*T/Q   x' = x / L
y
'=y/L   t' = tT / SL ^ 2
K = (K
'/B')L ^ 2 / T
epsilon = S / sigma
'where s=drawdown, T=transmissivity and S=storativity of the pumped aquifer. Also, t=time,
'K' = aquitard
hydraulic
conductivity, B'=aquitard saturated thickness, L=length chosen
'arbitrarily and sigma=aquitard specific yield or effective porosity at the location of the
'free surface. This uses the method of images, and nmax controls the number of image points
'used. (Increasing nmax by one adds four new image points to the image well system.) The
'prime superscript has been omitted for notational convenience.
Function
Eta_14(x, y, t, K, epsilon, alpha, beta, nmax)
If
t <= 0  # Then
Eta_14 = 0  #
Else
r = Sqr(x ^ 2 + y ^ 2)
Eta_14 = Eta_3(r, t, epsilon, K)
For
n = 0
To
nmax
x1 = alpha - (x + beta) + (2 * n + 1) * (alpha + beta)
x2 = alpha + (x + beta) + (2 * n + 1) * (alpha + beta)
x3 = beta - (alpha - x) + (2 * n + 1) * (alpha + beta)
x4 = beta + (alpha - x) + (2 * n + 1) * (alpha + beta)
r1 = Sqr(x1 ^ 2 + y ^ 2)
r2 = Sqr(x2 ^ 2 + y ^ 2)
r3 = Sqr(x3 ^ 2 + y ^ 2)
r4 = Sqr(x4 ^ 2 + y ^ 2)
Eta_14 = Eta_14 + Eta_3(r1, t, epsilon, K) + Eta_3(r2, t, epsilon, K) + _
Eta_3(r3, t, epsilon, K) + Eta_3(r4, t, epsilon, K)
Next
n
End
If
End
Function
"""

"""
below is from dispersion

'This module contains routines to calculate concentrations from contaminant sources in one,
'two and three-dimensional groundwater flow. Dimensionless independent variables are  defined
'as follows:
'
' (x',y',z',t',Dx',Dy',Dz',lambda')=(x/L,y/L,z/L,t*u/R*L,Dx/u*L,Dy/u*L,Dz/u*L,lambda*R*L/u)
'
'where u=pore velocity, L=length that may be chosen in any way, t=time, z=vertical coordinate,
'R=retardation factor, lambda=decay constant and u is in the positive x direction. Both x and
'y are horizontal coordinates, and dispersion coefficients are denoted by Dx, Dy and Dz.All
'input and output variables are dimensionless, and the prime superscript has been omitted for
'notational convenience. Dimensionless concentrations are defined at the beginning of each
'program.

'This calculates concentrations for a one-dimensional instantaneous source. The dimensionless
'concentration is c_1=c*sigma*L/M1 where c=concentration, sigma=porosity, L=length and M1=
'mass per unit area. (This area is normal to the direction of flow.)
Function c_1(x, t, Dx, lambda)
If t <= 0# Then
    c_1 = 0#
Else
    Pi = 3.141592654
    a = lambda * t + (x - t) ^ 2 / (4 * Dx * t)
    b = Sqr(4 * Pi * t * Dx)
    c_1 = Exp(-a) / b
End If
End Function

'This calculates concentrations for a two-dimensional instantaneous source. The dimensionless
'concentration is c_2=c*sigma*L^2/M2 where c=concentration, sigma=porosity, L=length and M2=
'mass per unit distance. [This distance is normal to the (x,y) plane.]
Function c_2(x, y, t, Dx, Dy, lambda)
If t <= 0# Then
    c_2 = 0#
Else
    Pi = 3.141592654
    a = lambda * t + (x - t) ^ 2 / (4 * Dx * t) + y ^ 2 / (4 * Dy * t)
    b = 4 * Pi * t * Sqr(Dx * Dy)
    c_2 = Exp(-a) / b
End If
End Function

'This calculates concentrations for a three-dimensional instantaneous source. The dimensionless
'concentration is c_3=c*sigma*L^3/M3 where c=concentration, sigma=porosity, L=length and M3=
'mass (in, for example, kg).
Function c_3(x, y, z, t, Dx, Dy, Dz, lambda)
If t <= 0# Then
    c_3 = 0#
Else
    Pi = 3.141592654
    a = lambda * t + (x - t) ^ 2 / (4 * Dx * t) + y ^ 2 / (4 * Dy * t) + z ^ 2 / (4 * Dz * t)
    b = 4 * Pi * t * Sqr(4 * Pi * t * Dx * Dy * Dz)
    c_3 = Exp(-a) / b
End If
End Function

'This calculates concentrations for a one-dimensional continuous source. The dimensionless
'concentration is c_4=c*sigma*u/R*dM1/dt where c=concentration, sigma=porosity, L=length,
'R=retardation factor and dM1/dt=constant mass injection rate of M1 at x=0. (M1=mass per unit
'area.This area is normal to the direction of flow.) This solution holds for
'-infinity < x < infinity.
Function c_4(x, t, Dx, lambda)
If t <= 0# Then
    c_4 = 0#
Else
    a = Abs(x) / Sqr(Dx)
    b = Sqr(lambda + 1 / (4 * Dx))
    c = x / (2 * Dx) - (a ^ 2 / (4 * t) + t * b ^ 2)
    d = a / (2 * Sqr(t)) + b * Sqr(t)
    e = x / (2 * Dx) - a * b
    f = a / (2 * Sqr(t)) - b * Sqr(t)
    c_4 = Exp(e) * Erfc(f) - Exp(c) * ExpErfc(d)
    c_4 = c_4 / (4 * b * Sqr(Dx))
End If
End Function

'This calculates concentrations for a two-dimensional continuous source. The dimensionless
'concentration is c_6=c*sigma*u*L/R*dM2/dt where c=concentration, sigma=porosity, L=length
'R=retardation factor and dM2/dt=constant mass injection rate of M2 at x=y=0. [M2=mass per unit
'distance normal to the (x,y) plane.]
Function c_6(x, y, t, Dx, Dy, lambda)
If t <= 0# Then
    c_6 = 0#
Else
    Pi = 3.141592654
    a = Sqr(x * x / Dx + y * y / Dy)
    b = Sqr(lambda + 1 / (4 * Dx))
    c_6 = Exp(x / (2 * Dx) - a * b) * Expw(a * a / (4 * t), a * b)
    c_6 = c_6 / (4 * Pi * Sqr(Dx * Dy))
End If
End Function

'This calculates concentrations for a three-dimensional continuous source. The dimensionless
'concentration is c_7=c*sigma*u*L^2/R*dM3/dt where c=concentration, sigma=porosity, L=length
'R=retardation factor and dM3/dt=constant mass injection rate of M3 at x=y=0. (M2=contaminate
'in kg.)
Function c_7(x, y, z, t, Dx, Dy, Dz, lambda)
If t <= 0# Then
    c_7 = 0#
Else
    Pi = 3.141592654
    a = Sqr(x * x / Dx + y * y / Dy + z * z / Dz)
    b = Sqr(lambda + 1 / (4 * Dx))
    c = a / (2 * Sqr(t)) + b * Sqr(t)
    d = a / (2 * Sqr(t)) - b * Sqr(t)
    e = x / (2 * Dx)
    c_7 = Exp(e - a * b) * Erfc(d) + Exp(e - a * a / (4 * t) - b * b * t) * ExpErfc(c)
    c_7 = c_7 / (8 * Pi * a * Sqr(Dx * Dy * Dz))
End If
End Function


'This calculates concentrations for a one-dimensional continuous source, where 0 < x < infinity
'and c=constant=c0 at x=0. Dimensionless variables are defined differently, as follows:
'
'        (c',x',t',lambda')=(c/c0,x*u/Dx,t*u^2/R*D,lambda*R*Dx/u^2)
'
'where dimensional variables are as defined previously. All input and output variables are
'dimensionless.
Function c_5(x, t, lambda)
If t <= 0# Then
    c_5 = 0#
Else
    a = x / (2 * Sqr(t))
    b = Sqr(t * (lambda + 1 / 4))
    c = a - b
    d = a + b
    c_5 = Exp(-c * c) * ExpErfc(d) + Erfc(c)
    c_5 = c_5 * Exp(x / 2 - 2 * a * b) / 2
End If
End Function

'This calculates concentrations for an instantaneous line source with finite length
'in one dimension. The dimensionless concentration is c_1_Line=c*sigma*L^2/M1 where
'c=concentration, sigma=porosity, L=arbitrarily chosen length and M1=total mass released
'instantaneously along the line (per unit area in the x-y plane). The line
'source extends from x_0 to x_M and is divided into M segments of equal
'length with one instantaneous, equal-strength source placed at the midpoint of each
'segment. Additional dimensionless variables follow:
'           x_0'=x_0/L        x_M'=x_M/L
'
Function c_1_Line(x, t, Dx, lambda, x_0, x_M, M)
If t <= 0# Then
    c_1_Line = 0#
Else
    c_1_Line = 0#
    For j = 1 To M
        xj = x_0 + (x_M - x_0) * (j - 1 / 2) / M
        c_1_Line = c_1_Line + c_1(x - xj, t, Dx, lambda)
    Next j
    c_1_Line = c_1_Line / M
End If
End Function

'This calculates concentrations for an instantaneous line source with finite length
'in two dimensions. The dimensionless concentration is c_2_Line=c*sigma*L^2/M2 where
'c=concentration, sigma=porosity, L=arbitrarily chosen length and M2=total mass released
'instantaneously along the line (per unit distance normal to the x-y plane). The line
'source extends from (x_0,y_0) to (x_M,y_M) and is divided into M segments of equal
'length with one instantaneous, equal-strength source placed at the midpoint of each
'segment. Additional dimensionless variables follow:
'   x_0'=x_0/L    y_0'=y_0/L     x_M'=x_M/L    y_M'=y_m/L
'
Function c_2_Line(x, y, t, Dx, Dy, lambda, x_0, y_0, x_M, y_M, M)
If t <= 0# Then
    c_2_Line = 0#
Else
    c_2_Line = 0#
    For j = 1 To M
        xj = x_0 + (x_M - x_0) * (j - 1 / 2) / M
        yj = y_0 + (y_M - y_0) * (j - 1 / 2) / M
        c_2_Line = c_2_Line + c_2(x - xj, y - yj, t, Dx, Dy, lambda)
    Next j
    c_2_Line = c_2_Line / M
End If
End Function

'This calculates concentrations for an instantaneous line source with finite length
'in three dimensions. The dimensionless concentration is c_3_Line=c*sigma*L^3/M3 where
'c=concentration, sigma=porosity, L=arbitrarily chosen length and M3=total mass released
'instantaneously along the line. The line source extends from (x_0,y_0,z_0) to
'(x_M,y_M,z_M) and is divided into M segments of equal length with one instantaneous,
'equal-strength source placed at the midpoint of each segment. Additional dimensionless
'variables follow:
'   x_0'=x_0/L    y_0'=y_0/L     z_0'=z_0/L   x_M'=x_M/L    y_M'=y_M/L    z_M'=z_M/L
'
Function c_3_Line(x, y, z, t, Dx, Dy, Dz, lambda, x_0, y_0, z_0, x_M, y_M, z_M, M)
If t <= 0# Then
    c_3_Line = 0#
Else
    c_3_Line = 0#
    For j = 1 To M
        xj = x_0 + (x_M - x_0) * (j - 1 / 2) / M
        yj = y_0 + (y_M - y_0) * (j - 1 / 2) / M
        zj = z_0 + (z_M - z_0) * (j - 1 / 2) / M
        c_3_Line = c_3_Line + c_3(x - xj, y - yj, z - zj, t, Dx, Dy, Dz, lambda)
    Next j
    c_3_Line = c_3_Line / M
End If
End Function

'This calculates concentrations for a continuous line source with finite length
'in one dimension. The dimensionless concentration is c_4_Line=c*sigma*u/R*dM1/dt where
'c=concentration, sigma=porosity, L=length, R=retardation factor and dM1/dt=constant mass
'injection rate of M1 from the line. (M1=mass per unit area.This area is normal to the
'direction of flow.) The line source extends from x_0 to x_M and is divided into M
'segments of equal length with one instantaneous, equal-strength source placed at the
'midpoint of each segment. Additional dimensionless variables follow:
'           x_0'=x_0/L        x_M'=x_M/L
'
Function c_4_Line(x, t, Dx, lambda, x_0, x_M, M)
If t <= 0# Then
    c_4_Line = 0#
Else
    c_4_Line = 0#
    For j = 1 To M
        xj = x_0 + (x_M - x_0) * (j - 1 / 2) / M
        c_4_Line = c_4_Line + c_4(x - xj, t, Dx, lambda)
    Next j
    c_4_Line = c_4_Line / M
End If
End Function

'This calculates concentrations for a continuous line source with finite length
'in two dimensions. The dimensionless concentration is c_6_Line=c*sigma*u*L/R*dM2/dt where
'c=concentration, sigma=porosity, L=length, R=retardation factor and dM2/dt=constant mass
'injection rate of M2 from the line. [M2=mass per unit distance normal to the (x,y) plane.]
'The line source extends from (x_0,y_0) to (x_M,y_M) and is divided into M
'segments of equal length with one instantaneous, equal-strength source placed at the
'midpoint of each segment. Additional dimensionless variables follow:
'      x_0'=x_0/L     y_0'=y_0/L    x_M'=x_M/L     y_M'=y_M/L
'
Function c_6_Line(x, y, t, Dx, Dy, lambda, x_0, y_0, x_M, y_M, M)
If t <= 0# Then
    c_6_Line = 0#
Else
    c_6_Line = 0#
    For j = 1 To M
        xj = x_0 + (x_M - x_0) * (j - 1 / 2) / M
        yj = y_0 + (y_M - y_0) * (j - 1 / 2) / M
        c_6_Line = c_6_Line + c_6(x - xj, y - yj, t, Dx, Dy, lambda)
    Next j
    c_6_Line = c_6_Line / M
End If
End Function

'This calculates concentrations for a continuous line source with finite length
'in three dimensions. The dimensionless concentration is c_7_Line=c*sigma*u*L^2/R*dM3/dt where
'c=concentration, sigma=porosity, L=length, R=retardation factor and dM2/dt=constant mass
'injection rate of M3 from the line. The line source extends from (x_0,y_0,z_0) to
'(x_M,y_M,z_M) and is divided into M segments of equal length with one instantaneous,
'equal-strength source placed at the midpoint of each segment. Additional dimensionless
'variables follow:
'   x_0'=x_0/L    y_0'=y_0/L   z_0'=z_0/L   x_M'=x_M/L   y_M'=y_M/L    z_M'=z_M/L
'
Function c_7_Line(x, y, z, t, Dx, Dy, Dz, lambda, x_0, y_0, z_0, x_M, y_M, z_M, M)
If t <= 0# Then
    c_7_Line = 0#
Else
    c_7_Line = 0#
    For j = 1 To M
        xj = x_0 + (x_M - x_0) * (j - 1 / 2) / M
        yj = y_0 + (y_M - y_0) * (j - 1 / 2) / M
        zj = z_0 + (z_M - z_0) * (j - 1 / 2) / M
        c_7_Line = c_7_Line + c_7(x - xj, y - yj, z - zj, t, Dx, Dy, Dz, lambda)
    Next j
    c_7_Line = c_7_Line / M
End If
End Function

'This calculates concentrations for an instantaneous three-dimensional source in an
'aquifer with a saturated depth L. Therefore, L cannot be chosen arbitrarily, and the
'dimensionless concentration is c_8=c*sigma*L^3/M3 where M3=contaminant mass. The
'contaminant source has the dimensionless coordinates (x',y',z')=(0,0,z_0/L) at t=0.
Function c_8(x, y, z, t, Dx, Dy, Dz, lambda, z_0)
If t <= 0 Then
    c_8 = 0
Else
    Pi = 3.141592654
    c_8 = 1
    n = 0
    Do
        n = n + 1
        alpha = n * Pi
        er = 2 * Exp(-t * Dz * alpha ^ 2)
        term = er * Cos(alpha * z_0) * Cos(alpha * z)
        c_8 = c_8 + term
    Loop Until er < 0.000001
    c_8 = c_8 * c_2(x, y, t, Dx, Dy, lambda)
End If
End Function

'This calculates concentrations for an instantaneous two-dimensional source in an
'aquifer with a saturated depth L. Therefore, L cannot be chosen arbitrarily, and the
'dimensionless concentration is c_8=c*sigma*L^2/M2 where M2=contaminant mass per unit
'distance in the y direction. The contaminant source has the dimensionless coordinates
'(x',z')=(0,z_0/L) at t=0.
Function c_9(x, z, t, Dx, Dz, lambda, z_0)
If t <= 0 Then
    c_9 = 0
Else
    Pi = 3.141592654
    c_9 = 1
    n = 0
    Do
        n = n + 1
        alpha = n * Pi
        er = 2 * Exp(-t * Dz * alpha ^ 2)
        term = er * Cos(alpha * z_0) * Cos(alpha * z)
        c_9 = c_9 + term
    Loop Until er < 0.000001
    c_9 = c_9 * c_1(x, t, Dx, lambda)
End If
End Function

'This calculates concentrations for a continuous three-dimensional source in an
'aquifer with a saturated depth L. Therefore, L cannot be chosen arbitrarily, and the
'dimensionless concentration is c_10=c*sigma*u*L^2/(R*dM3/dt) where dM3/dt=contaminant
'mass injection rate. The contaminant source has the dimensionless coordinates
'(x',y',z')=(0,0,z_0/L).
Function c_10(x, y, z, t, Dx, Dy, Dz, lambda, z_0)
If t <= 0 Then
    c_10 = 0
Else
    Pi = 3.141592654
    c_10 = c_6(x, y, t, Dx, Dy, lambda)
    n = 0
    Do
        n = n + 1
        alpha = n * Pi
        er = 2 * c_6(x, y, t, Dx, Dy, lambda + Dz * alpha ^ 2)
        term = er * Cos(alpha * z_0) * Cos(alpha * z)
        c_10 = c_10 + term
    Loop Until er < 0.000001
End If
End Function

'This calculates concentrations for a continuous two-dimensional source in an
'aquifer with a saturated depth L. Therefore, L cannot be chosen arbitrarily, and the
'dimensionless concentration is c_11=c*sigma*u*L/(R*dM2/dt) where dM2/dt=contaminant
'mass injection rate. The contaminant source has the dimensionless coordinates
'(x',z')=(0,z_0/L).
Function c_11(x, z, t, Dx, Dz, lambda, z_0)
If t <= 0 Then
    c_11 = 0
Else
    Pi = 3.141592654
    c_11 = c_4(x, t, Dx, lambda)
    n = 0
    Do
        n = n + 1
        alpha = n * Pi
        er = 2 * c_4(x, t, Dx, lambda + Dz * alpha ^ 2)
        term = er * Cos(alpha * z_0) * Cos(alpha * z)
        c_11 = c_11 + term
    Loop Until er < 0.000001
End If
End Function

'This calculates concentrations for an instantaneous line source with finite length
'in three dimensions. The dimensionless concentration is c_8_Line=c*sigma*L^3/M3 where
'c=concentration, sigma=porosity, L=aquifer saturated depth and M3=total mass released
'instantaneously along the line. The line source extends from (x_0,y_0,z_0) to
'(x_M,y_M,z_M) and is divided into M segments of equal length with one instantaneous,
'equal-strength source placed at the midpoint of each segment. Additional dimensionless
'variables follow:
'   x_0'=x_0/L    y_0'=y_0/L     z_0'=z_0/L   x_M'=x_M/L    y_M'=y_M/L    z_M'=z_M/L
'
Function c_8_Line(x, y, z, t, Dx, Dy, Dz, lambda, x_0, y_0, z_0, x_M, y_M, z_M, M)
If t <= 0# Then
    c_8_Line = 0#
Else
    c_8_Line = 0#
    For j = 1 To M
        xj = x_0 + (x_M - x_0) * (j - 1 / 2) / M
        yj = y_0 + (y_M - y_0) * (j - 1 / 2) / M
        zj = z_0 + (z_M - z_0) * (j - 1 / 2) / M
        c_8_Line = c_8_Line + c_8(x - xj, y - yj, z, t, Dx, Dy, Dz, lambda, zj)
    Next j
    c_8_Line = c_8_Line / M
End If
End Function

'This calculates concentrations for an instantaneous line source with finite length
'in two dimensions. The dimensionless concentration is c_9_Line=c*sigma*L^2/M2 where
'c=concentration, sigma=porosity, L=aquifer saturated depth and M2=total mass per unit
'distance in the y direction released instantaneously along a line in the (x,z) plane.
'The line source extends from (x_0,z_0) to (x_M,z_M) and is divided into M segments
'of equal length with one instantaneous,equal-strength source placed at the midpoint of each
'segment. Additional dimensionless variables follow:
'        x_0'=x_0/L      z_0'=z_0/L   x_M'=x_M/L      z_M'=z_M/L
'
Function c_9_Line(x, z, t, Dx, Dz, lambda, x_0, z_0, x_M, z_M, M)
If t <= 0# Then
    c_9_Line = 0#
Else
    c_9_Line = 0#
    For j = 1 To M
        xj = x_0 + (x_M - x_0) * (j - 1 / 2) / M
        yj = y_0 + (y_M - y_0) * (j - 1 / 2) / M
        zj = z_0 + (z_M - z_0) * (j - 1 / 2) / M
        c_9_Line = c_9_Line + c_9(x - xj, z, t, Dx, Dz, lambda, zj)
    Next j
    c_9_Line = c_9_Line / M
End If
End Function

'This calculates concentrations for a continuous line source with finite length
'in three dimensions. The dimensionless concentration is c_10_Line=c*sigma*u*L^2/(R*dM3/dt)
'where c=concentration, sigma=porosity, L=aquifer saturated depth and dM3/dt=total contaminant
'mass flow rate injected continuously along the line. The line source extends from
'(x_0,y_0,z_0) to (x_M,y_M,z_M) and is divided into M segments of equal length with one
'continuous, equal-strength source placed at the midpoint of each segment. Additional
'dimensionless variables follow:
'   x_0'=x_0/L    y_0'=y_0/L     z_0'=z_0/L   x_M'=x_M/L    y_M'=y_M/L    z_M'=z_M/L
'
Function c_10_Line(x, y, z, t, Dx, Dy, Dz, lambda, x_0, y_0, z_0, x_M, y_M, z_M, M)
If t <= 0# Then
    c_10_Line = 0#
Else
    c_10_Line = 0#
    For j = 1 To M
        xj = x_0 + (x_M - x_0) * (j - 1 / 2) / M
        yj = y_0 + (y_M - y_0) * (j - 1 / 2) / M
        zj = z_0 + (z_M - z_0) * (j - 1 / 2) / M
        c_10_Line = c_10_Line + c_10(x - xj, y - yj, z, t, Dx, Dy, Dz, lambda, zj)
    Next j
    c_10_Line = c_10_Line / M
End If
End Function

'This calculates concentrations for a continuous line source with finite length
'in two dimensions. The dimensionless concentration is c_11_Line=c*sigma*u*L/(R*dM2/dt
'where c=concentration, sigma=porosity, L=aquifer saturated depth and dM2/dt=contaminant
'mass flow rate injected continuously along the line. The line source extends from
'(x_0,z_0) to (x_M,z_M) and is divided into M segments of equal length with one continuous,
'equal-strength source placed at the midpoint of each segment. Additional dimensionless
'variables follow:
'         x_0'=x_0/L     z_0'=z_0/L   x_M'=x_M/L     z_M'=z_M/L
'
Function c_11_Line(x, z, t, Dx, Dz, lambda, x_0, z_0, x_M, z_M, M)
If t <= 0# Then
    c_11_Line = 0#
Else
    c_11_Line = 0#
    For j = 1 To M
        xj = x_0 + (x_M - x_0) * (j - 1 / 2) / M
        yj = y_0 + (y_M - y_0) * (j - 1 / 2) / M
        zj = z_0 + (z_M - z_0) * (j - 1 / 2) / M
        c_11_Line = c_11_Line + c_11(x - xj, z, t, Dx, Dz, lambda, zj)
    Next j
    c_11_Line = c_11_Line / M
End If
End Function

"""

"""
below is from math functions 

 'This module contains  routines to calculate mathematical functions.


'To compute the gamma function
Function gamma(x)
    gamma = Exp(GAMMLN(x))
End Function
                
'To compute the modified Bessel function I0(x) for 0<x<infinity.
Function BessI0(x)
A0 = 1
a1 = 3.5156229
a2 = 3.0899424
A3 = 1.2067492
A4 = 0.2659732
A5 = 0.0360768
A6 = 0.0045813
B0 = 0.39894228
B1 = 0.01328592
B2 = 0.00225319
B3 = -0.00157565
B4 = 0.00916281
B5 = -0.02057706
B6 = 0.02635537
B7 = -0.01647633
B8 = 0.00392377
If x <= 3.75 Then
    t = (x / 3.75) ^ 2
    BessI0 = A0 + t * (a1 + t * (a2 + t * (A3 + t * (A4 + t * (A5 + t * A6)))))
Else
    t = 3.75 / x
    BessI0 = B0 + t * (B1 + t * (B2 + t * (B3 + t * (B4 + t * (B5 + t * (B6 + t * (B7 + t * B8)))))))
    BessI0 = BessI0 * Exp(x) / Sqr(x)
End If
End Function

'To compute the modified Bessel function I0(x) multiplied by exp(-x) for 0<x<infinity.
Function ExpBessI0(x)
A0 = 1
a1 = 3.5156229
a2 = 3.0899424
A3 = 1.2067492
A4 = 0.2659732
A5 = 0.0360768
A6 = 0.0045813
B0 = 0.39894228
B1 = 0.01328592
B2 = 0.00225319
B3 = -0.00157565
B4 = 0.00916281
B5 = -0.02057706
B6 = 0.02635537
B7 = -0.01647633
B8 = 0.00392377
If x <= 3.75 Then
    t = (x / 3.75) ^ 2
    ExpBessI0 = Exp(-x) * (A0 + t * (a1 + t * (a2 + t * (A3 + t * (A4 + t * (A5 + t * A6))))))
Else
    t = 3.75 / x
    ExpBessI0 = B0 + t * (B1 + t * (B2 + t * (B3 + t * (B4 + t * (B5 + t * (B6 + t * (B7 + t * B8)))))))
    ExpBessI0 = ExpBessI0 / Sqr(x)
End If
End Function
   
   
'To compute the modified Bessel function I1(x) for 0<x<infinity.
Function BessI1(x)
A0 = 0.5
a1 = 0.87890594
a2 = 0.51498869
A3 = 0.15084934
A4 = 0.02658733
A5 = 0.00301532
A6 = 0.00032411
B0 = 0.39894228
B1 = -0.03988024
B2 = -0.00362018
B3 = 0.00163801
B4 = -0.01031555
B5 = 0.02282967
B6 = -0.02895312
B7 = 0.01787654
B8 = -0.00420059
If x <= 3.75 Then
    t = (x / 3.75) ^ 2
    BessI1 = A0 + t * (a1 + t * (a2 + t * (A3 + t * (A4 + t * (A5 + t * A6)))))
    BessI1 = BessI1 * x
Else
    t = 3.75 / x
    BessI1 = B0 + t * (B1 + t * (B2 + t * (B3 + t * (B4 + t * (B5 + t * (B6 + t * (B7 + t * B8)))))))
    BessI1 = BessI1 * Exp(x) / Sqr(x)
End If
End Function

'To compute the modified Bessel function I1(x) multiplied by exp(-x) for 0<x<infinity.
Function ExpBessI1(x)
A0 = 0.5
a1 = 0.87890594
a2 = 0.51498869
A3 = 0.15084934
A4 = 0.02658733
A5 = 0.00301532
A6 = 0.00032411
B0 = 0.39894228
B1 = -0.03988024
B2 = -0.00362018
B3 = 0.00163801
B4 = -0.01031555
B5 = 0.02282967
B6 = -0.02895312
B7 = 0.01787654
B8 = -0.00420059
If x <= 3.75 Then
    t = (x / 3.75) ^ 2
    ExpBessI1 = A0 + t * (a1 + t * (a2 + t * (A3 + t * (A4 + t * (A5 + t * A6)))))
    ExpBessI1 = ExpBessI1 * x * Exp(-x)
Else
    t = 3.75 / x
    ExpBessI1 = B0 + t * (B1 + t * (B2 + t * (B3 + t * (B4 + t * (B5 + t * (B6 + t * (B7 + t * B8)))))))
    ExpBessI1 = ExpBessI1 / Sqr(x)
End If
End Function
                
'To calculate the modified Bessel function K0(x) for 0<x<infinity.
Function BessK0(x)
A0 = -0.57721566
a1 = 0.4227842
a2 = 0.23069756
A3 = 0.0348859
A4 = 0.00262698
A5 = 0.0001075
A6 = 0.0000074
B0 = 1.25331414
B1 = -0.07832358
B2 = 0.02189568
B3 = -0.01062446
B4 = 0.00587872
B5 = -0.0025154
B6 = 0.00053208
If x <= 2 Then
    t = (x / 2) ^ 2
    BessK0 = A0 + t * (a1 + t * (a2 + t * (A3 + t * (A4 + t * (A5 + t * A6)))))
    BessK0 = BessK0 - Application.Ln(x / 2) * BessI0(x)
Else
    t = 2 / x
    BessK0 = B0 + t * (B1 + t * (B2 + t * (B3 + t * (B4 + t * (B5 + t * B6)))))
    BessK0 = BessK0 * Exp(-x) / Sqr(x)
End If
End Function
                    
                    
'To compute the modified Bessel function K1(x) for 0<x<infinity.
Function BessK1(x)
A0 = 1
a1 = 0.15443144
a2 = -0.67278579
A3 = -0.18156897
A4 = -0.01919402
A5 = -0.00110404
A6 = -0.00004686
B0 = 1.25331414
B1 = 0.23498619
B2 = -0.0365562
B3 = 0.01504268
B4 = -0.00780353
B5 = 0.00325614
B6 = -0.00068245
If x <= 2 Then
    t = (x / 2) ^ 2
    BessK1 = A0 + t * (a1 + t * (a2 + t * (A3 + t * (A4 + t * (A5 + t * A6)))))
    BessK1 = BessK1 / x + Application.Ln(x / 2) * BessI1(x)
Else
    t = 2 / x
    BessK1 = B0 + t * (B1 + t * (B2 + t * (B3 + t * (B4 + t * (B5 + t * B6)))))
    BessK1 = BessK1 * Exp(-x) / Sqr(x)
End If
End Function


        
        
'To compute the Bessel function J(n,x) for n=integer and 0<x<infinity.
Function BessJ(n, x)
a = Sqr((Abs(x) + 7) ^ 2 / 2 + n ^ 2 / 4)
L = Int(2 + 2 * Int(a))
BJ = 1
d = 1
r = 1
Do
    r = x / (2 * L - r * x)
    If L <= n Then
        BJ = BJ * r
    End If
    Frac = L / 2 - Int(L / 2)
    d = d * r + 2 * Frac
    L = L - 1
Loop Until L = 0
BessJ = BJ / (2 * d - 1)
End Function
         

'To compute the Bessel function Y0(x)for 0<x<infinity.
Function BessY0(x)
A0 = 0.36746691
a1 = 0.60559366
a2 = -0.74350384
A3 = 0.25300117
A4 = -0.04261214
A5 = 0.00427916
A6 = -0.00024846
B0 = 0.79788456
B1 = -0.00000077
B2 = -0.0055274
B3 = -0.00009512
B4 = 0.00137237
B5 = -0.00072805
B6 = 0.00014476
c0 = -0.78539816
c1 = -0.04166397
c2 = -0.00003954
c3 = 0.00262573
C4 = -0.00054125
C5 = -0.00029333
C6 = 0.00013558
If x <= 3 Then
    t = (x / 3) ^ 2
    BessY0 = A0 + t * (a1 + t * (a2 + t * (A3 + t * (A4 + t * (A5 + t * A6)))))
    BessY0 = BessY0 + (2 / Application.Pi()) * Application.Ln(x / 2) * BessJ(0, x)
Else
    t = 3 / x
    P0 = B0 + t * (B1 + t * (B2 + t * (B3 + t * (B4 + t * (B5 + t * B6)))))
    p1 = x + c0 + t * (c1 + t * (c2 + t * (c3 + t * (C4 + t * (C5 + t * C6)))))
    BessY0 = P0 * Sin(p1) / Sqr(x)
End If
End Function
               
               
'To compute the Bessel function BessY1(x)for 0<x<infinity.
Function BessY1(x)
A0 = -0.6366198
a1 = 0.2212091
a2 = 2.1682709
A3 = -1.3164827
A4 = 0.3123951
A5 = -0.0400976
A6 = 0.0027873
B0 = 0.79788456
B1 = 0.00000156
B2 = 0.01659667
B3 = 0.00017105
B4 = -0.00249511
B5 = 0.00113653
B6 = -0.00020033
c0 = -2.35619449
c1 = 0.12499612
c2 = 0.0000565
c3 = -0.00637879
C4 = 0.00074348
C5 = 0.00079824
C6 = -0.00029166
If x <= 3 Then
    t = (x / 3) ^ 2
    BessY1 = A0 + t * (a1 + t * (a2 + t * (A3 + t * (A4 + t * (A5 + t * A6)))))
    BessY1 = BessY1 / x + (2 / Application.Pi()) * Application.Ln(x / 2) * BessJ(1, x)
Else
    t = 3 / x
    P0 = B0 + t * (B1 + t * (B2 + t * (B3 + t * (B4 + t * (B5 + t * B6)))))
                       
    p1 = x + c0 + t * (c1 + t * (c2 + t * (c3 + t * (C4 + t * (C5 + t * C6)))))
    BessY1 = P0 * Sin(p1) / Sqr(x)
End If
End Function
                       

'To compute the exponential integral Exp1(x) for 0<x<infinity.
Function Exp1(x)
A0 = -0.57721566
a1 = 0.99999193
a2 = -0.24991055
A3 = 0.05519968
A4 = -0.00976004
A5 = 0.00107857
B0 = 0.2677737343
B1 = 8.6347608925
B2 = 18.059016973
B3 = 8.5733287401
B4 = 1
c0 = 3.9584969228
c1 = 21.0996530827
c2 = 25.6329561486
c3 = 9.5733223454
C4 = 1
If x <= 1 Then
    Exp1 = -Log(x) + A0 + x * (a1 + x * (a2 + x * (A3 + x * (A4 + x * A5))))
Else
    p1 = B0 + x * (B1 + x * (B2 + x * (B3 + x * B4)))
    P2 = c0 + x * (c1 + x * (c2 + x * (c3 + x * C4)))
    Exp1 = (p1 / P2) * Exp(-x) / x
End If
End Function
        
       
'To compute the complimentary error function for 0<=Abs(x)<infinity.
Function Erfc1(x)
    u = Abs(x)
    p = 0.3275911
    a1 = 0.254829592
    a2 = -0.284496736
    A3 = 1.421413741
    A4 = -1.453152027
    A5 = 1.061405429
    t = 1 / (1 + p * u)
    a = t * (a1 + t * (a2 + t * (A3 + t * (A4 + t * A5)))) * Exp(-u * u)
    If x >= 0 Then
        Erfc1 = a
    Else
        Erfc1 = 2 - a
    End If
End Function

'To compute the error function for 0<=Abs(x)<infinity.
Function Erf(x)
    Erf = 1 - Erfc1(x)
End Function

'To compute the product exp(x^2)*Erfc1(x).
Function ExpErfc(x)
    If x >= 5 Then
        Pi = 3.14159265
        ExpErfc = (1 - 0.5 / x ^ 2 + 0.75 / x ^ 4) / (x * Sqr(Pi))
    Else
        p = 0.3275911
        a1 = 0.254829592
        a2 = -0.284496736
        A3 = 1.421413741
        A4 = -1.453152027
        A5 = 1.061405429
        t = 1 / (1 + p * x)
        ExpErfc = t * (a1 + t * (a2 + t * (A3 + t * (A4 + t * A5))))
    End If
End Function
    
'To compute the exponential integral of order n for n >=0.
Function ExpInt(n, x)
If n = 0 Then
    ExpInt = Exp(-x) / x
ElseIf n = 1 Then
    ExpInt = Exp1(x)
ElseIf (n > 1) And (x <= 5) Then
    a = Exp1(x)
    For i = 2 To n
        a = (Exp(-x) - x * a) / (i - 1)
    Next i
    ExpInt = a
ElseIf (n > 1) And (x > 5) Then
    N1 = Int(x)
    t = x + N1
    a = 1 + N1 / t ^ 2 + N1 * (N1 - 2 * x) / t ^ 4 + N1 * (6 * x ^ 2 - 8 * N1 * x + N1 ^ 2) / t ^ 6
    a = a * Exp(-x) / t
    If n <= N1 Then
        i = N1
        Do While i > n
            i = i - 1
            a = (Exp(-x) - i * a) / x
        Loop
        ExpInt = a
    Else
        i = N1
        Do While i < n
            i = i + 1
            a = (Exp(-x) - x * a) / (i - 1)
        Loop
        ExpInt = a
    End If
End If
End Function

'To compute the Hantush leaky-aquifer function W(x,y). Append routines to
'   compute Exp1(x),ExpInt(n,x),BessK0(x),BessI0(x)and BessI1(x).
Function w(x, y)
If x = 0 Then
    w = 2 * BessK0(y)
Else
    r = 1
    t = y ^ 2 / (4 * x)
    b = 2 * x
    If y <= b Then
        w = 0
        n = 0
        Do
            term = r * ExpInt(n + 1, x)
            w = w + term
            n = n + 1
            r = r * (-t) / n
        Loop Until Abs(term) < 0.0000000001
    Else
        w = 2 * BessK0(y)
        n = 0
        Do
            term = r * ExpInt(n + 1, t)
            w = w - term
            n = n + 1
            r = r * (-x) / n
        Loop Until Abs(term) < 0.0000000001
    End If
End If
End Function

'To compute the Hantush leaky-aquifer function scaled with the exponential of its second
'argument, Exp(y)*w(x,y).
Function Expw(x, y)
If y < 5 Then
    Expw = Exp(y) * w(x, y)
Else
    Pi = 3.141592654
    g = Sqr(x) - y / (2 * Sqr(x))
    i = Erfc1(g) * Sqr(Pi) / 2
    S = i
    term = 1#
    n = 0
    Do
        n = n + 1
        i = i * (2 * n - 1) / 2 + Exp(-g ^ 2) * g ^ (2 * n - 1) / 2
        term = -term * (2 * n - 1) / (n * (4 * y))
        er = term * i
        S = S + er
    Loop Until Abs(er) < 0.0001
    Expw = S * Sqr(2 / y)
End If
End Function

'To compute the partial derivative of W(x,y) with respect to x.
Function Wx(x, y)
If x = 0# Then
    Wx = 0#
Else
    Wx = -Exp(-x - y ^ 2 / (4 * x)) / x
End If
End Function

'To compute the partial derivative of W(x,y) with respect to y.
Function Wy(x, y)
If x = 0# Then
    x = 0.000001
End If
If y = 0# Then
    y = 0.000001
End If
r = 1
t = y ^ 2 / (4 * x)
b = 2 * x
If y <= b Then
    Wy = 0#
    n = 0
    Do
        term = r * ExpInt(n + 2, x)
        Wy = Wy + term
        n = n + 1
        r = r * (-t) / n
    Loop Until Abs(term) < 0.0000000001
    Wy = -Wy * y / (2 * x)
Else
    Wy = 0#
    n = 0
    Do
        term = r * ExpInt(n, t)
        Wy = Wy + term
        n = n + 1
        r = r * (-x) / n
    Loop Until Abs(term) < 0.0000000001
    Wy = -2 * BessK1(y) + Wy * y / (2 * x)
End If
End Function

'To compute the Airy function Ai(x) for both positive and negative arguments.
Function Ai(x)
    j = 1
    Q = 1
    If 5 < Abs(x) Then GoTo 2
    t = 0.614926628
    u = x * Sqr(3) / (2 * t * Application.Pi())
    a = t - u
    z = x ^ 3 / 3
1:
    t = t * z / (j * (3 * j - 1))
    u = u * z / (j * (3 * j + 1))
    a = a + t - u
    j = j + 1
    If j < (2 + 3 * Abs(x)) Then GoTo 1
    a = a / Sqr(3)
    GoTo 6
2:
    If x > 0 Then GoTo 3
    Q = -1
    x = -x
3:
    a = 1 / Sqr(2 * Application.Pi() * Sqr(x))
    t = a
    b = a
    u = a
    z = 2 * x ^ 1.5 / 3
4:
    b = b * (j - 1 + 1 / (7.2 * j)) / (2 * z)
    a = -a * (j - 1 + 1 / (7.2 * j)) / (2 * z)
    t = t + a
    u = u + b
    j = j + 1
    If j < (3 + 90 / z) Then GoTo 4
    If Q < 0 Then GoTo 5
    a = t * Exp(-z) / Sqr(2)
    GoTo 6
5:
    a = t * Cos(z) + u * Sin(z)
6:
    Ai = a
End Function


' Calculates the value of the natural log of the gamma function using
' an algotithm given on page 157 of "Numerical Recipes'
    Function GAMMLN(XX)
    Dim COF(6)
    COF(1) = 76.18009173
    COF(2) = -86.50532033
    COF(3) = 24.01409822
    COF(4) = -1.231739516
    COF(5) = 0.00120858003
    COF(6) = -0.00000536382
    STP = 2.50662827465
    HALF = 0.5
    ONE = 1#
    FPF = 5.5
    x = XX - ONE
    TMP = x + FPF
    TMP = (x + HALF) * Log(TMP) - TMP
    SER = ONE
    j = 1
    Do
        x = x + ONE
        SER = SER + COF(j) / x
        j = j + 1
    Loop Until j = 7
    GAMMLN = TMP + Log(STP * SER)
    End Function
    
' This calculates the beta function. The file GAMMLN(XX) must be attached.
    Function beta(z, w)
    beta = Exp(GAMMLN(z) + GAMMLN(w) - GAMMLN(z + w))
    End Function
    
   
' This calculates the incomplete beta function Ix(a,b) using an
' algorithm given on page 167 of "Numerical Recipes". The files
' GAMMLN(XX) and BETACF(A,B,X) must be attached. Restrictions on the
' arguments are 0<=X<=1, 0<A and 0<B.
    Function BETAI(a, b, x)
    If x = 0# Then
        BR = 0#
    Else
        If x = 1# Then
            BR = 0#
        Else
            BT = Exp(GAMMLN(a + b) - GAMMLN(a) - GAMMLN(b) + a * Log(x) + b * Log(1# - x))
        End If
    End If
    If x < (a + 1#) / (a + b + 2#) Then
        BETAI = BT * BETACF(a, b, x) / a
    Else
        BETAI = 1# - BT * BETACF(b, a, 1# - x) / b
    End If
    End Function

' This is used by BETAI(A,B,X) to calculate a continued fraction
' for the incomplete beta function.
    Function BETACF(a, b, x)
    EPS = 0.0000003
    AM = 1#
    BM = 1#
    AZ = 1#
    QAB = a + b
    QAP = a + 1#
    QAM = a - 1#
    BZ = 1# - QAB * x / QAP
    M = 1
    Do
    EM = M
    TEM = EM + EM
    d = EM * (b - M) * x / ((QAM + TEM) * (a + TEM))
    AP = AZ + d * AM
    BP = BZ + d * BM
    d = -(a + EM) * (QAB + EM) * x / ((a + TEM) * (QAP + TEM))
    APP = AP + d * AZ
    BPP = BP + d * BZ
    AOLD = AZ
    AM = AP / BPP
    BM = BP / BPP
    AZ = APP / BPP
    BZ = 1#
    Loop Until Abs(AZ - AOLD) < EPS * Abs(AZ)
    BETACF = AZ
    End Function
    
' This calculates the incomplete gamma function P(a,x) by using
' an algorithm given on page 161 of "Numerical Recipes". The file
' GAMMLN(XX) must be attached. Restrictions on the arguments
' are X>=0 and A>0.
    Function GAMMP(a, x)
    EPS = 0.0000003
    GLN = GAMMLN(a)
    If x = 0# Then
        GAMMP = 0#
    ElseIf x < (a + 1) Then
        AP = a
        sum1 = 1# / a
        DEL = sum1
        Do
        AP = AP + 1
        DEL = DEL * x / AP
        sum1 = sum1 + DEL
        Loop Until Abs(DEL) < Abs(sum1) * EPS
        GAMMP = sum1 * Exp(-x + a * Log(x) - GLN)
    Else
        GOLD = 0#
        A0 = 1#
        a1 = x
        B0 = 0#
        B1 = 1#
        FAC = 1#
        N1 = 1
        Do
        GOLD1 = GOLD
        AN = N1
        ANA = AN - a
        A0 = (a1 + A0 * ANA) * FAC
        B0 = (B1 + B0 * ANA) * FAC
        ANF = AN * FAC
        a1 = x * A0 + ANF * a1
        B1 = x * B0 + ANF * B1
        FAC = 1# / a1
        g = B1 * FAC
        GOLD = g
        N1 = N1 + 1
        Loop Until Abs(g - GOLD1) / g < EPS
        GAMMP = 1# - Exp(-x + a * Log(x) - GLN) * g
    End If
    End Function


' This calculates the incomplete gamma function Q(a,x) by using
' an algorithm given on page 162 of "Numerical RPecipes". The file
' GAMMLN.for must be attached.Restrictions on the arguments
' are X>=0 and A>0.
    Function GAMMQ(a, x)
    EPS = 0.0000003
    GLN = GAMMLN(a)
    If x = 0# Then
        GAMMQ = 1#
    ElseIf x < (a + 1) Then
        AP = a
        sum1 = 1# / a
        DEL = sum1
        Do
        AP = AP + 1
        DEL = DEL * x / AP
        sum1 = sum1 + DEL
        Loop Until Abs(DEL) < Abs(sum1) * EPS
        GAMMQ = 1 - sum1 * Exp(-x + a * Log(x) - GLN)
    Else
        GOLD = 0#
        A0 = 1#
        a1 = x
        B0 = 0#
        B1 = 1#
        FAC = 1#
        N1 = 1
        Do
        GOLD1 = GOLD
        AN = N1
        ANA = AN - a
        A0 = (a1 + A0 * ANA) * FAC
        B0 = (B1 + B0 * ANA) * FAC
        ANF = AN * FAC
        a1 = x * A0 + ANF * a1
        B1 = x * B0 + ANF * B1
        FAC = 1# / a1
        g = B1 * FAC
        GOLD = g
        N1 = N1 + 1
        Loop Until Abs(g - GOLD1) / g < EPS
        GAMMQ = Exp(-x + a * Log(x) - GLN) * g
    End If
    End Function
    

'This computes the incomplete gamma function, gamma*, using an algorithm given
'on page 441 of Spanier and Oldham.
Function Gammastar(nu, x)
    g = 1
    p = 1
    nu = 1 + nu
    Do While nu < 2
        g = g * x
        p = p * nu + g
        nu = 1 + nu
    Loop
    j = Int(5 * (3 + Abs(x)) / 2)
    f = 1 / (j + nu - x)
    Do While j > 0
        j = j - 1
        f = (f * x + 1) / (j + nu)
    Loop
    p = p + f * g * x
    a = 2 / (7 * nu ^ 2)
    b = 2 / (3 * nu ^ 2)
    g = (1 - a * (1 - b)) / (30 * nu ^ 2)
    g = (g - 1) / (12 * nu) - nu * (Log(nu) - 1)
    Gammastar = p * Exp(g - x) * Sqr(nu / (2 * Application.Pi()))
End Function

'This calculates the integral of funct(x) from x=a to  x=b using the trapezoidal
'rule with n slices. funct(x) must be computed with another user-defined program.
Function trapz(a, b, n)
    delta = (b - a) / n
    trapz = 0#
    y_1 = funct(a)
    For i = 1 To n
        y_2 = funct(a + i * delta)
        trapz = trapz + delta * (y_1 + y_2) / 2
        y_1 = y_2
    Next i
End Function

'This calculates the integral of funct(x) from x=a to  x=b using Simpson's one-third
'rule with n slices. n must be an even number, and funct(x) must be computed with
'another user-defined program.
Function Simp(a, b, n)
    n = Application.Even(n)
    delta = (b - a) / n
    Simp = 0#
    y_1 = funct(a)
    For i = 1 To n - 1 Step 2
        y_2 = funct(a + i * delta)
        y_3 = funct(a + (i + 1) * delta)
        Simp = Simp + delta * (y_1 + 4 * y_2 + y_3) / 3
        y_1 = y_3
    Next i
End Function

'This computes funct(x)
Function funct(x)
    funct = Exp(-x)
End Function

'This uses the Stehfest algorithm to invert a Laplace transform. See Stehfest,H.
'(1970) "Numerical inversion of Laplace transforms," Comm. ACM, Vol.13, No.1,
'pages 47-49 and 624. The approximation has the nature of an asymptotic series,
'in which n is the number of terms used. Thus, n shuld not exceed about 20, and
'n = 10 is probably sufficient for most applications. Note that n must be an even
'integer and that t = time. Beware: this works well for some transforms but very
'poorly for others. Furthermore, a value for n that is slightly too large can
'cause a dramatic decrease in accuracy.
Function Lapinv(t, n)
n = Application.Even(n)
Lapinv = 0#
For i = 1 To n
    Lapinv = Lapinv + Stehcoef(i, n) * transform(i * Log(2) / t)
Next i
Lapinv = Lapinv * Log(2) / t
End Function

'This calculates the coefficient c(i) in the Stehfest algorithm for Laplace transform
'inversion. The integer n must be even.
Function Stehcoef(i, n)
M = Application.Round(n / 2, 0)
upperlimit = Application.Min(i, M)
lowerlimit = Application.RoundDown((i + 1) / 2, 0)
Stehcoef = 0#
For K = lowerlimit To upperlimit
    num = Application.Fact(2 * K) * K ^ M
    denom = Application.Fact(M - K) * Application.Fact(K) * Application.Fact(K - 1) * _
        Application.Fact(i - K) * Application.Fact(2 * K - i)
    Stehcoef = Stehcoef + num / denom
Next K
Stehcoef = Stehcoef * (-1) ^ (i + M)
End Function

'This calculates a numerical value for the Laplace transform of a function for use in
'Lapinv(t,n). p = transform parameter.
Function transform(p)
transform = 2 * BessK0(2 * Sqr(p)) / p
End Function

'To compute ecosh(x)=exp(-|x|)*cosh(x)=(1+exp(-2|x|))/2 for -infinity<x<infinity. In words,
' this scales cosh(x) with exp(-|x|) for all real values of x.
Function ecosh(x)
ecosh = (1 + Exp(-2 * Abs(x))) / 2
End Function

'To compute esinh(x)=exp(-|x|)*sinh(x)=sign(x)*(1-exp(-2|x|)/2 for -infinity<x<infinity.
'In words, this scales sinh(x) with exp(-|x|) for all real values of x.
Function esinh(x)
esinh = Sgn(x) * (1 - Exp(-2 * Abs(x))) / 2
End Function

'This calculates the integral of funct(x)*cos(c*x) from x=a to  x=b using a specialized
'quadrature formula. This routine is designed to work for -infinity<c<infinity and reduces
'to Simpson's rule when c=0. The routine uses n slices, where n must be an even number.
'The function funct(x) must be computed with another user-defined program.
Function IntegralCos(a, b, c, n)
n = Application.Even(n)
delta = (b - a) / n
xi = c * delta
If xi <> 0 Then
    c0 = (Sin(xi) - xi * Cos(xi)) * 4 / xi ^ 3
    c1 = (xi * Cos(xi) - Sin(xi)) / xi ^ 2
    c2 = (2 * xi * Cos(xi) + (xi ^ 2 - 2) * Sin(xi)) / xi ^ 3
Else
    c0 = 4 / 3
    c1 = 0#
    c2 = 1 / 3
End If
IntegralCos = 0#
y_1 = funct(a)
For i = 1 To n - 1 Step 2
    cx = c * (a + i * delta)
    coef_1 = c2 * Cos(cx) - c1 * Sin(cx)
    coef_2 = c0 * Cos(cx)
    coef_3 = c1 * Sin(cx) + c2 * Cos(cx)
    y_2 = funct(a + i * delta)
    y_3 = funct(a + (i + 1) * delta)
    IntegralCos = IntegralCos + delta * (coef_1 * y_1 + coef_2 * y_2 + coef_3 * y_3)
    y_1 = y_3
Next i
End Function

'This calculates the integral of funct(x)*sin(c*x) from x=a to  x=b using a specialized
'quadrature formula. This routine is designed to work for -infinity<c<infinity and reduces
'to zero when c=0. The routine uses n slices, where n must be an even number. The function
'funct(x) must be computed with another user-defined program.
Function IntegralSin(a, b, c, n)
n = Application.Even(n)
delta = (b - a) / n
xi = c * delta
If xi <> 0 Then
    c0 = (Sin(xi) - xi * Cos(xi)) * 4 / xi ^ 3
    c1 = (xi * Cos(xi) - Sin(xi)) / xi ^ 2
    c2 = (2 * xi * Cos(xi) + (xi ^ 2 - 2) * Sin(xi)) / xi ^ 3
Else
    c0 = 0#
    c1 = 0#
    c2 = 0#
End If
IntegralSin = 0#
y_1 = funct(a)
For i = 1 To n - 1 Step 2
    cx = c * (a + i * delta)
    coef_1 = c2 * Sin(cx) + c1 * Cos(cx)
    coef_2 = c0 * Sin(cx)
    coef_3 = -c1 * Cos(cx) + c2 * Sin(cx)
    y_2 = funct(a + i * delta)
    y_3 = funct(a + (i + 1) * delta)
    IntegralSin = IntegralSin + delta * (coef_1 * y_1 + coef_2 * y_2 + coef_3 * y_3)
    y_1 = y_3
Next i
End Function


"""

"""
below is from scale dependencies

'This module contains routines to calculate concentrations for instantaneous and steady-state
'sources when dispersivities increase directly with the first power of distance downstream.
'All input and output variables are dimensionless. Independent variables for instantaneous
'sources are defined as follows:
'  x'=xR/(ut); y'=yR/(ut); z'=zR/(ut)
'where R=retardation coefficient and u=pore velocity.
'Independent variables for steady-state sources are defined as follows:
'                     y'=y/x; z'=z/x
'Dispersion coefficients are given by D1=Epsilon1*u*x, D2=Epsilon2*u*x and D3=Epsilon3*u*x.
'The prime superscripts are omitted for notational convenience. These results are taken from
'Hunt, B. (1998) "Contaminant source solutions with scale-dependent dispersivities", Jnl.
'Hydrologic Engrg., ASCE, Vol.3, No.4, 268-275. Note that Eq. (62) and Fig. 6 contain errors
'in this paper.




'This computes concentrations for an instantaneous source in one dimension. The dimensionless
'concentration is given by
'           C'=C*Sigma*u*t/M1
Function Cinstx_1(x, Epsilon1)
If x = 0# Then
    Cinstx_1 = 0#
Else
    e = Epsilon1
    x = x / e
    a = (1 + e / 12 + e ^ 2 / 288 - 139 * e ^ 3 / 51840 - 571 * e ^ 4 / 2488320) * Sqr(2 * e * Application.Pi())
    Cinstx_1 = Exp(-x + 1 / e + Log(x * e) / e) / a
End If
End Function


'This computes concentrations for an instantaneous source in two dimensions. The dimensionless
'concentration is given by
'           C'=C*Sigma*(u*t)^2/M2
Function Cinstx_2(x, y, Epsilon1, Epsilon2)
    x = x / Epsilon1
    y = y / Sqr(Epsilon1 * Epsilon2)
    n = 200
    delta = 20 / n
    integral = 0#
    p = 0#
    For i = 1 To n Step 2
        term = F_2(x, y, Epsilon1, p) + 4 * F_2(x, y, Epsilon1, p + delta) + F_2(x, y, Epsilon1, p + 2 * delta)
        integral = integral + term
        p = p + 2 * delta
    Next i
    integral = integral * delta / 3
    Cinstx_2 = integral * x ^ (1 / Epsilon1) / (Application.Pi() * Sqr(Epsilon1 * Epsilon2) * gamma(1 / Epsilon1))
End Function


'The function to be integrated.
Function F_2(x, y, Epsilon1, p)
    If p = 0 Then
        a = 1
        b = Exp(-x)
    Else
        a = (Application.Sinh(p) / p) ^ (1 + 1 / Epsilon1)
        b = Exp(-x * p / Application.Tanh(p)) * Cos(p * y)
    End If
    F_2 = b / a
End Function


'This computes concentrations for an instantaneous source in three dimensions. The dimensionless
'concentration is given by
'           C'=C*Sigma*(u*t)^3/M3
Function Cinstx_3(x, y, z, Epsilon1, Epsilon2, Epsilon3)
x = x / Epsilon1
y = y / Sqr(Epsilon1 * Epsilon2)
z = z / Sqr(Epsilon1 * Epsilon3)
n = 200
    delta = 20 / n
    integral = 0#
    p = 0#
    For i = 1 To n Step 2
        term = f_3(x, y, z, Epsilon1, p) + 4 * f_3(x, y, z, Epsilon1, p + delta) + f_3(x, y, z, Epsilon1, p + 2 * delta)
        integral = integral + term
        p = p + 2 * delta
    Next i
    integral = integral * delta / 3
    Cinstx_3 = integral * x ^ (1 / Epsilon1) / (2 * Application.Pi() * gamma(1 / Epsilon1) * Epsilon1 * Sqr(Epsilon2 * Epsilon3))
End Function


'The function to be integrated.
Function f_3(x, y, z, Epsilon1, p)
    If p = 0 Then
        a = 1
        b = 0#
    Else
        a = (Application.Sinh(p) / p) ^ (1 + 1 / Epsilon1)
        b = p * BessJ(0, p * Sqr(y ^ 2 + z ^ 2)) * Exp(-x * p / Application.Tanh(p))
    End If
    f_3 = b / a
End Function


'This computes concentrations for a steady-state source in two dimensions. The dimensionless
'concentration is given by
'           C'=C*Sigma*u*x/(R*dM2/dt)
Function Cstdyx_2(y, Epsilon1, Epsilon2)
    y = y * Sqr(Epsilon1 / Epsilon2)
    a = gamma((1 + 1 / Epsilon1) / 2) * Sqr(Epsilon1 / Epsilon2)
    a = a / (gamma(1 / 2) * gamma(1 / (2 * Epsilon1)))
    r = (1 + y ^ 2) ^ ((1 + 1 / Epsilon1) / 2)
    Cstdyx_2 = a / r
End Function


'This computes concentrations for a steady-state source in three dimensions. The dimensionless
'concentration is given by
'           C'=C*Sigma*u*x^2/(R*dM3/dt)
Function Cstdyx_3(y, z, Epsilon1, Epsilon2, Epsilon3)
    y = y * Sqr(Epsilon1 / Epsilon2)
    z = z * Sqr(Epsilon1 / Epsilon3)
    a = 2 * Application.Pi() * Sqr(Epsilon2 * Epsilon3)
    r = (1 + y ^ 2 + z ^ 2) ^ (1 + 1 / (2 * Epsilon1))
    Cstdyx_3 = 1 / (a * r)
End Function


'This computes concentrations for the boundary condition c(0,t)=c0 for 0<t<infinity. The
'initial condition is c(x,0)=0 for 0<x<infinity, and the dispersion coefficient is epsilon*u*x.
'All input and output variables are dimensionless and are defined as follows:
'             c'=c/c0; t'=t*u/(x*R); Lambda'=R*Lambda*x/u
'where R=retardation factor and Lambda=radioactive decay constant.
Function c_cont(epsilon, t, lambda)
    e = epsilon
    L = lambda
If t = 0# Then
    c_cont = 0#
Else
    integral = 0#
    n = Int(5 * t / Sqr(e) + 1)
    n = 2 * n
    delta = t / n
    x = 0#
    a = 1 + e / 12 + e ^ 2 / 288 - 139 * e ^ 3 / 51840 - 571 * e ^ 4 / 2488320
    a = a * Sqr(2 * e * Application.Pi())
    For i = 1 To n Step 2
        If x = 0# Then
            y1 = 0#
        Else
            y1 = Exp(-L * x - (1 / x + Log(x) - 1) / e) / x
        End If
        y2 = Exp(-L * (x + delta) - (1 / (x + delta) + Log(x + delta) - 1) / e) / (x + delta)
        y3 = Exp(-L * (x + 2 * delta) - (1 / (x + 2 * delta) + Log(x + 2 * delta) - 1) / e) / (x + 2 * delta)
        term = y1 + 4 * y2 + y3
        integral = integral + term
        x = x + 2 * delta
    Next i
    integral = integral * delta / 3
    c_cont = integral / a
End If
End Function


'This computes concentration for a pulse. The dispersion coefficient is epsilon*u*x,
'and the flow is one-dimensional with the boundary condition c(0,t)=c0 for 0<t<t0 and
'c(0,t)=0 for t0<t<infinity. The initial condition is c(x,0)=0 for 0<x<infinity.
'All input and output variables are dimensionless and are defined as follows:
'              c'=c/c0; t'=t*u/(x*R); t_0'=t_0*u/(x*R); Lambda'=R*Lambda*x/u
'where t_0 = pulse duration time, R=retardation factor and Lambda=radioactive decay constant.
Function c_pulse(epsilon, t, t_0, lambda)
    If t <= t_0 Then
        c_pulse = c_cont(epsilon, t, lambda)
    Else
        c_pulse = c_cont(epsilon, t, lambda) - c_cont(epsilon, t - t_0, lambda)
    End If
End Function


'This uses a perturbation approximation for Epsilon1<<1 and Epsilon2<<1 to calculate an unsteady
'solution for contaminant transport from a pit of length 2*L. The concentration at the pit is c0.
'Definitions for the dimensionless variables follow:
'             c'=c/c0; x'=x/L; y'=y/L; t'=tu/(RL); Lambda'=Lambda*R*L/u
'The longitudinal and lateral dispersion coefficients are Epsilon1*u*x and Epsilon2*u*x,
'respectively. All input and output variables are dimensionless, and all prime superscripts have
'been omitted for notational convenience. The following reference discusses this solution:
'Hunt, B. 2002. "Scale-Dependent Dispersion from a Pit", Jnl. Hydrologic Engrg.,
'Amer. Soc. Civ. Engrs., Vol. 7, No. 2, pp. 168-174.
Function Pit_unstdy2(x, y, t, Epsilon1, Epsilon2, lambda)
    If x = 0 Then
        If Abs(y) <= 1 Then
            Pit_unstdy2 = 1
        Else
            Pit_unstdy2 = 0
        End If
    Else
            a = x * Sqr(2 * Epsilon2)
            b = t * Sqr(2 * Epsilon1)
            c = t * Sqr(2 * Epsilon2)
        If x < t Then
            d = Exp(-x * lambda) * (Erf((1 - y) / a) + Erf((1 + y) / a)) / 2
            d = d - Exp(-t * lambda) * Erfc((t - x) / b) * (Erf((1 - y) / c) + Erf((1 + y) / c)) / 4
            Pit_unstdy2 = d
        Else
            d = Exp(-t * lambda) * Erfc((x - t) / b) * (Erf((1 - y) / c) + Erf((1 + y) / c)) / 4
            Pit_unstdy2 = d
        End If
    End If
End Function

"""