import numpy as np

def cos(a):
    return np.cos(np.radians(a))

def sin(a):
    return np.sin(np.radians(a))

def tan(a):
    return np.tan(np.radians(a))

def acos(a):
    return np.degrees(np.arccos(a))

def asin(a):
    return np.degrees(np.arcsin(a))

def atan(a):
    return np.degrees(np.arctan(a))

def atan2(a, b):
    return np.degrees(np.arctan2(a, b))