% -*- coding: utf-8 -*-
% # Adaptive Thiele Continued Fractions #

% This Jupyter notebook contains Maple code (licensed by the MIT license below) to explore the use of the Adaptive Thiele Continued Fraction method of rational approximation to approximate the Lambert $W$ function near the branch point using 
% \begin{equation}
% f(s) = \left\{\begin{array}{cc}
% W_0\left(-\exp(-1)(1-s^2/2)\right) & s \ge 0
% \\
% W_{-1}\left(-\exp(-1)(1-s^2/2)\right) & s < 0  
% \end{array}\right.
% \end{equation}

% © Robert M. Corless 2026 (MIT Release licence)
%
% Permission is hereby granted, free of charge, to any person obtaining a copy of this software and associated documentation files (the "Software"), to deal in the Software without restriction, including without limitation the rights to use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies of the Software, and to permit persons to whom the Software is furnished to do so, subject to the following conditions:
%
% The above copyright notice and this permission notice shall be included in all copies or substantial portions of the Software.
%
% THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
%

%
% Script to approximate Lambert W near the branch point
%
clf
n = 117;
w0s = linspace(-1,0,n);
x0 = w0s.*exp(w0s);
s0 = sqrt(2*(1+exp(1.)*x0));
w1s = linspace(-8,-1,n);
x1 = w1s.*exp(w1s);
s1 = -sqrt(2*(1+exp(1.)*x1));
% ws(1) is the same as w1s(n), both give s=0 and w=-1
% x0(1) is the same as x1(n);
s = [s1(1:n-1),0,s0(2:n)];
w = [w1s(1:n-1),-1,w0s(2:n)];
x = w.*exp(w);
[r,pol,res,zer,aa, zz, err, ssc] = tcf( w, s, 'tol', 2.0e-13, 'eqtm', 1 );
err(end-3:end)
ss = linspace(-sqrt(2),sqrt(2), 10240*n);
wss = r(ss);
xss = -exp(-1.)*(1-ss.^2/2);
figure(1), plot( s, w, 'k.', ss, wss, 'b')
berr = xss - wss.*exp(wss);
figure(2), semilogy( ss, abs(berr), 'r.')
figure(3), plot(real(pol),imag(pol),'k*')
[mm, ~] = size(err);
figure(4), semilogy( 1:mm, err, 'kx')
grid
figure(5), plot( ss, berr,'k.')
axis([-sqrt(2),sqrt(2),-3.e-14,3.e-14])
ax = gca
ax.FontSize = 24
xlabel('s', 'FontSize',24)
ylabel('error','FontSize',24)
grid


