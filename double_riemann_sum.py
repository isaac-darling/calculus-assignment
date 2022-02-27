#Isaac Darling, February 2021
#Write a function to simplify the process for computing double riemann sums

#returns a numeric value based on mathematical function f()
def f(x, y):
    return 100-x**2+6*x*y-3*y**2

"""
params:
   func=mathematical function to estimate integral for,
   bounds=list of two lists of length 2 containg coords for bottom left and top right,
   m=number of x axis segments,
   n=number of y axis segments
returns midpoint double riemann sum, but the ranges can be altered to do left or right sums
    simply replace list comprehension with range(m) or range(n) for left
    and range(1,m+1) or range(1,n+1) for right
"""
def drs(func, bounds, m, n):
    summation = 0
    dx = (bounds[1][0] - bounds[0][0])/m
    dy = (bounds[1][1] - bounds[0][1])/n
    for i in [x+0.5 for x in range(m)]:
        for j in [x+0.5 for x in range(n)]:
            summation+=func(i*dx, j*dy)
    return summation*dx*dy

print(drs(f, [[0, 0], [1, 1]], 3, 3))