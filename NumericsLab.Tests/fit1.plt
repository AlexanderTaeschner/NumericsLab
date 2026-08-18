u(x) = x
v(x) = 16 - x
min(a,b) = a < b ? a : b
w(x) = min(u(x),v(x))
f(x) = x1 + u(x)/(v(x)*x2 + w(x)*x3)

x1 = 1
x2 = 1
x3 = 1


set fit noerrorscaling errorvariables covariancevariables
fit f(x) 'fit1.dat' via x1,x2,x3

plot 'fit1.dat', f(x)