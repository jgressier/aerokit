import numpy as np
import inspect

@np.vectorize
def secant_solve(f, y, x, tol=0.00001):
    """
    secant_solve(f, y, x, tol=0.00001)
      secant method (approximate newton) searching for x with f(x)=y
    """
    itmax = 20
    it    = 1
    old_x = x
    old_y = f(old_x)
    new_x = x*0.999 
    #while (it < itmax) and np.any(abs(new_x-old_x) > np.maximum(tol,tol*abs(new_x))):
    while (it < itmax) and (abs(new_x-old_x) > max(tol,tol*abs(new_x))):
        new_y = f(new_x)
        next_x = new_x + (new_x-old_x)/(new_y-old_y)*(y-new_y)
        old_x = new_x
        old_y = new_y
        new_x = next_x
        it += 1
    if it == itmax:
        print("!!! max number of iterations (%2i) reached"%(it))
        callseq = ''
        for i in range(len(inspect.stack())):
            callseq += '/'+inspect.stack()[i][3]
        print('    from',callseq)
    return new_x
