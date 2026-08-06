import numpy as np
from aerokit.common._ode import RK4, _rkf45

def test_RK4_basic_functionality():
    # Testing the RK4 function with a simple ODE dy/dx = y, whose solution is y = e^x
    F = lambda x, y: y  # The derivative function
    x0, y0 = 0, 1  # Initial condition: y(0) = 1
    xStop = 1  # Stop at x = 1
    h = 0.1  # Step size
    X, Y = RK4(F, x0, y0, xStop, h)
    # Ensure outputs are numpy arrays
    assert isinstance(X, np.ndarray), "X should be a numpy array"
    assert isinstance(Y, np.ndarray), "Y should be a numpy array"
    # Check if the length of X and Y matches the expected number of steps
    assert len(X) == len(Y), "X, and Y array lengths mismatch"
    # Verify approximate solution at x = 1 (expecting around e^1 ~ 2.718)
    assert np.isclose(Y[-1], y0*np.exp(xStop), atol=1e-8), "Final value mismatch in RK4"

# def test_RK4_edge_case_zero_step():
#     # Testing RK4 with a zero step size
#     F = lambda x, y: y
#     x0, y0 = 0, 1
#     xStop = 1
#     h = 0  # Zero step size

#     X, Y = RK4(F, x0, y0, xStop, h)
#     # With zero step, should only have initial values
#     assert len(X) == 1 and len(Y) == 1, "Expected single initial point in output"
#     assert X[0] == x0 and Y[0] == y0, "Initial values do not match in zero step test"

def test_rkf45_basic_functionality():
    # Testing the RKF45 function with the same ODE dy/dx = y
    F = lambda x, y: np.array([y])
    x0, y0 = 0, np.array([1])  # Initial condition in array form
    h = 0.1  # Initial step size
    # Call _rkf45 method directly with these parameters
    K, e = _rkf45(F, x0, y0, h)
    # Check if K is a numpy array and has the expected shape
    assert isinstance(K, np.ndarray), "K should be a numpy array"
    assert K.shape == (1,), "Unexpected shape of output K in rkf45"
    assert np.isclose(K, y0*(np.exp(h)-1.), atol=1.e-12), "bad one step integration"
    assert e < 1.e-8, "too large estimation of error"

def test_rkf45_correctness_on_simple_ode():
    # Verify RKF45's correctness on dy/dx = y, over one step
    F = lambda x, y: np.array([y])
    x0, y0 = 0, np.array([1])  # y(0) = 1
    h = 0.1

    K = _rkf45(F, x0, y0, h)
    # Verify expected increment from initial condition using a close tolerance
    assert np.isclose(K[0], h * y0, atol=0.1).all(), "Mismatch in RKF45 increment for simple ODE"
