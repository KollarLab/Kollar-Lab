from functools import wraps

def freeze(cls):
    cls.__frozen = False

    def frozensetattr(self, key, value):
        if self.__frozen and not hasattr(self, key):
            raise AttributeError("Class {} is frozen. Cannot set {} = {}"
                  .format(cls.__name__, key, value))
        else:
            object.__setattr__(self, key, value)

    def init_decorator(func):
        @wraps(func)
        def wrapper(self, *args, **kwargs):
            func(self, *args, **kwargs)
            self.__frozen = True
        return wrapper

    cls.__setattr__ = frozensetattr
    cls.__init__ = init_decorator(cls.__init__)

    return cls






###############
#shamelessly stolen and modified elipse fitting functions
#################
import numpy as np
from numpy.linalg import eig, inv
    
def _fitEllipse(x,y):
    '''This is the raw version that we stole from online, but it has parameter wierdness
    and branch cut issues, so I'm fixing it to return normal stuff and keep major and minor
    axes in a fixed order'''
    x = x[:,np.newaxis]
    y = y[:,np.newaxis]
    D =  np.hstack((x*x, x*y, y*y, x, y, np.ones_like(x)))
    S = np.dot(D.T,D)
    C = np.zeros([6,6])
    C[0,2] = C[2,0] = 2; C[1,1] = -1
    E, V =  eig(np.dot(inv(S), C))
    n = np.argmax(np.abs(E))
    a = V[:,n]
    return a    
    
def ellipse_center(a):
    b,c,d,f,g,a = a[1]/2, a[2], a[3]/2, a[4]/2, a[5], a[0]
    num = b*b-a*c
    x0=(c*d-b*f)/num
    y0=(a*f-b*d)/num
    return np.array([x0,y0])



def ellipse_axis_length( a ):
    b,c,d,f,g,a = a[1]/2, a[2], a[3]/2, a[4]/2, a[5], a[0]
    up = 2*(a*f*f+c*d*d+g*b*b-2*b*d*f-a*c*g)
    down1=(b*b-a*c)*( (c-a)*np.sqrt(1+4*b*b/((a-c)*(a-c)))-(c+a))
    down2=(b*b-a*c)*( (a-c)*np.sqrt(1+4*b*b/((a-c)*(a-c)))-(c+a))
    res1=np.sqrt(up/down1)
    res2=np.sqrt(up/down2)
    return np.array([res1, res2])
    
#def ellipse_angle_of_rotation2( a ):
#    b,c,d,f,g,a = a[1]/2, a[2], a[3]/2, a[4]/2, a[5], a[0]
#    if b == 0:
#        if a > c:
#            return 0
#        else:
#            return np.pi/2
#    else:
#        if a > c:
#            return np.arctan(2*b/(a-c))/2
#        else:
#            return np.pi/2 + np.arctan(2*b/(a-c))/2   
#        
def ellipse_angle_of_rotation( a ):
    '''modified from ellipse_angle_of_rotation2 to handle some residual branch cut
    issues with inverse trig functions'''
    b,c,d,f,g,a = a[1]/2, a[2], a[3]/2, a[4]/2, a[5], a[0]
    if b == 0:
        if a > c:
            out = 0
        else:
            out =  np.pi/2
    else:
        if a > c:
            out =  np.arctan(2*b/(a-c))/2
        else:
            out =  np.pi/2 + np.arctan(2*b/(a-c))/2           
        
    return out

def make_elipse(axes, phi, center,  numPoints):
    a, b = axes
    thetas = np.linspace(0,2*np.pi, numPoints)
    
    xx = center[0] + a*np.cos(thetas)*np.cos(phi) - b*np.sin(thetas)*np.sin(phi)
    yy = center[1] + a*np.cos(thetas)*np.sin(phi) + b*np.sin(thetas)*np.cos(phi)
    return xx, yy        
        
def fitEllipse(Is, Qs, verbose = False):
    ''' fits an ellipse using least squares method from online
    
    Then processes the fit parameters into more sensble things
    fixing the order of major and minor axes and fixing branch cuts
    in inverse trig functions
    
    Returns:
    [major axis, minor axis] , [x0, y0], angle of major axis to x axis    
    
    angle in radians'''
    
    #stolen online, lest squares elipse fit
    a = _fitEllipse(Is,Qs)
    center = ellipse_center(a)
    phi = ellipse_angle_of_rotation(a)
    axes = ellipse_axis_length(a)
    
    
    #fix the angle of the stig because this fit function 
    # is having domain issues
    #This mixer seems to have major axis along pi/4, roughly
    if phi > np.pi/4:
        phi = phi - np.pi/2
        
    if axes[1] > axes[0]:
        #major axis is second
        axes = [axes[1], axes[0]]
        phi = phi + np.pi/2
        
    if phi < 0:
        phi = phi +np.pi
    
#    xOffset = center[0]
#    yOffset = center[1]
    stigAngle =  180*phi/np.pi #degrees
    stig = (axes[1]-axes[0])/np.mean(axes)
    
    if verbose:
        print("    ")
        print("center = ",  np.round(center,3))
        print("angle of rotation = " + str(np.round( stigAngle, 3)) + ' degrees')
        print("axes = ", np.round(axes,3))
        print("stig = ", np.round(stig,3))
    
    return axes, center, phi        
        
###############        
        
        
        
        
        
        
        