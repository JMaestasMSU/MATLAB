# Bisection Method: computes approximate solution of f(x)=0
# Input: function f; a,b such that f(a)*f(b)<0,
# and tolerance tol
# Output: Approximate solution x such that f(x)=0
def sign(x):
    if x < 0:
        return -1
    elif x > 0:
        return 1
    else:
        return 0


def mybisect(f, a, b, tol):
    fa = f(a)
    fb = f(b)
    # note that we do not compute f(a)*f(b) since that could cause an overflow
    if sign(fa) * sign(fb) >= 0:
        # Interval does not bracket a root; raise so caller can handle it
        raise ValueError('f(a)*f(b) < 0 not satisfied for interval [{}, {}]'.format(a, b))

    # Main bisection loop: stop when half-interval width <= tol
    while (b - a) / 2.0 > tol:  # this is half the width of the interval
        c = (a + b) / 2.0
        fc = f(c)  # this is our only function evaluation this iteration
        if fc == 0:  # c is an exact root (rare for floats)
            return c
        if sign(fc) * sign(fa) < 0:  # root is in [a,c]
            b = c
            # fb = fc  # not needed unless you later use fb
        else:  # root is in [c,b]
            a = c
            fa = fc

    # Return midpoint of final interval as the approximate root
    return (a + b) / 2.0


def f(x):
    # corrected polynomial to match MATLAB: x^4 - 2x^3 - 4x^2 + 4x + 4
    return x**4 - 2.0*x**3 - 4.0*x**2 + 4.0*x + 4.0

tol = 1.0e-2
try:
    x = mybisect(f, -2, -1, tol)
    print(x)
    print(f(x))
except Exception as e:
    print('Error:', e)