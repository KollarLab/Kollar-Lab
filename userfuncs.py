from functools import wraps
from scipy.signal import butter,filtfilt
import numpy as np
from numpy.linalg import eig, inv
import pickle
import os
from datetime import datetime
from datetime import date
from skimage.measure import EllipseModel

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

def make_elipse(axes, center, phi,  numPoints):
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
    ecc = np.sqrt(abs(axes[0]**2-axes[1]**2))/max(axes)
    
    if verbose:
        print("    ")
        print("center = ",  np.round(center,3))
        print("angle of rotation = " + str(np.round( stigAngle, 3)) + ' degrees')
        print("axes = ", np.round(axes,3))
        print("stig = ", np.round(stig,3))
        print("ecc = ", np.round(ecc,3))
    
    return axes, center, phi        

def fit_ell_martin(Is, Qs, verbose = False):
    
    a_points = np.array(list(zip(Is,Qs)))   
    ell = EllipseModel()
    ell.estimate(a_points)
    
    xc, yc, a, b, theta = ell.params
    axes = [a,b]
    center = [xc, yc]
    phi = theta*180/np.pi
    ecc = np.sqrt(abs(axes[0]**2-axes[1]**2))/max(axes)
    if verbose:
        print("center = ",  np.round(center,3))
        print("angle of rotation = ",  np.round(phi,3))
        print("axes = ", np.round(axes,3))
        print("ecc = ", np.round(ecc,3))
    return axes, center, phi, ecc

###############
    
##########################
#Filters
##########################
def butter_lowpass_filter(data, cutoff, fs, order):
    normal_cutoff = cutoff / (fs/2)
    # Get the filter coefficients 
    b, a = butter(order, normal_cutoff, btype='low', analog=False)
    y = filtfilt(b, a, data)
    return y
      
#############
#figures
############       

def savefig(fig, name, path = '', png = False):
    '''
    save an interactive pickle version of a figure
    fig = figure object
    name = name of file, not including extension
    path = folder where the figure should end up
    png = boolean for if you want a png copy too
    
    '''
    
    if png:
        saveName = name + '.png'
        pathStr= os.path.join(path, saveName) 
        fig.savefig(pathStr, dpi = 200, transparent  = False)
        return
    
    if name[-4:] == '.pkl':
        saveName = name
    else:
        saveName = name + '.pkl'
    pathStr= os.path.join(path, saveName) 
    pickle.dump(fig, open(pathStr, 'wb'))

    return

def loadfig(path):
    figx = pickle.load(open(path, 'rb'))
    figx.show()
    return figx

def SaveInst(instruments):
    HWsettings = {}
    for inst in instruments.keys():
        HWsettings[inst] = instruments[inst].settings
    return HWsettings

def makedict(vars, localdict):
    return {var:localdict[var] for var in vars}

def SaveFull(path, name, variables, localdict, expsettings={}, instruments={}, figures=[], saveHWsettings=True):
    
    if name[-4:] == '.pkl':
        saveName = name
    else:
        saveName = name + '.pkl'
        
    pathStr               = os.path.join(path,saveName)

    toSave                = {}
    toSave['Data']        = makedict(variables, localdict)
    toSave['ExpSettings'] = expsettings
    if saveHWsettings:
        toSave['HWSettings']  = SaveInst(instruments)
    toSave['Figures']     = figures

    pickle.dump(toSave, open(pathStr, 'wb'))

def LoadFull(path):
    fullData    = pickle.load(open(path,'rb'))

    figures     = fullData['Figures']
    expsettings = fullData['ExpSettings']
    hwsettings  = fullData['HWSettings']
    data        = fullData['Data']

    return [data, expsettings, hwsettings, figures]

def timestamp():
    date  = datetime.now()
    stamp = date.strftime('%Y%m%d_%H%M%S')
    return stamp

def saveDir(settings):
    exp_globals  = settings['exp_globals']
    exp_settings = settings['exp_settings']
    
    today = date.today().strftime('%Y%m%d')
    
    root    = exp_globals['root_folder']
    project = exp_globals['project_name']
    device  = exp_globals['device_name']
    
    try:
        meas_type = exp_settings['meas_type']
    except:
        meas_type = exp_settings['spec']['meas_type']
    
    fullpath = os.path.join(root, project, device, meas_type, today)

    try:
        os.makedirs(fullpath)
    except:
        print('Dir already exists')
    return fullpath

#def saveDir(project_dir, meas_type):
#    today = date.today().strftime('%Y%m%d')
#    fullpath = os.path.join(project_dir, meas_type, today)
#    try:
#        os.makedirs(fullpath)
#    except:
#        print('Dir already exists')
#    return fullpath

def reset_local_vars(local_dict, global_dict, vars_to_save):
    
    save_list = vars_to_save
    
    for name in list(global_dict):
        if name in save_list:
            continue
        else:
            try:
                del global_dict[name]
            except:
                pass
            try:
                del local_dict[name]
            except:
                pass