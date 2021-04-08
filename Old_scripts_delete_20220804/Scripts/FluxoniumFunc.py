### creating code to functionalize basic fluxoonium functions

import numpy as np
import scipy.sparse as ssp
import scipy.sparse.linalg
import matplotlib.pyplot as plt
from matplotlib.widgets import Slider, Button, RadioButtons
from matplotlib.widgets import CheckButtons
import scipy.io
import matplotlib.gridspec as gridspec
import h5py
from mat4py import loadmat
from mpl_toolkits.mplot3d import Axes3D

###############################################################################
###############################################################################
###############################################################################
##########

#plt.rcParams['lines.linewidth'] = 2
#plt.rcParams['lines.markersize'] = 4

##########

# General constants

q_e = 1.602e-19  # electron charge in C
h = 6.626e-34  # Planck constant in Js
hbar = h / (2 * np.pi)  # Planck constant / 2p in Js
phi_0 = hbar / (2 * q_e)  # flux quantum / 2p in Wb

# beta_phi  = 0.2 # original given value for Beta

beta_phi = 0.20  # hopeful value
f_res = 7.8

g_phi = 0.125 * beta_phi * f_res #g_PHI = g_factor*charge_matrix element in GHz

tol = 1e-8  # diagonalization tolerance
keig = 8  # number of required ekets
N_grid = 251  # number of points in each direction (odd)
max_grid = 4 * np.pi  # maximum range of each coordinate
grid_pts = np.linspace(-max_grid, max_grid, N_grid, endpoint=True,
                       dtype='float64')
grid_x = keig
grid_y = 1
legend = []
scaling = 15


###############################################################################
###############################################################################
### Defining functions for finding the basic fluxonium parameters
###############################################################################

#from numba import jit
#@jit(nopython=True)
#

def write2vtk(matrix, fname='Tree_flux.vtk'):
    (sx, sy, sz) = matrix.shape
    norm = np.linalg.norm(matrix)
    lines = \
        '''# vtk DataFile Version 2.0
Volume example
ASCII
DATASET STRUCTURED_POINTS
DIMENSIONS %d %d %d
ASPECT_RATIO 1 1 1
ORIGIN 0 0 0
POINT_DATA %d
SCALARS matlab_scalars float 1
LOOKUP_TABLE default
''' \
        % (sx, sy, sz, matrix.size)
    with open(fname, 'w') as f:
        f.write(lines)
        for ix in range(sz):
            v = np.ravel(matrix[:, :, ix], order='f')
            v = ['%1.5f' % x for x in 100 * v / norm]
            line = ' '.join(v)
            f.write(line + '\n')


def hessian(x):
    """
    Calculate the hessian matrix with finite differences
    Parameters:
       - x : ndarray
    Returns:
       an array of shape (x.dim, x.ndim) + x.shape
       where the array[i, j, ...] corresponds to the second derivative x_ij
    """

    # Taken from https://stackoverflow.com/questions/
    # 31206443/numpy-second-derivative-of-a-ndimensional-array

    x_grad = np.gradient(x)
    hessian = np.empty((x.ndim, x.ndim) + x.shape, dtype=x.dtype)
    for (k, grad_k) in enumerate(x_grad):

        # iterate over dimensions
        # apply gradient again to every component of the first derivative.

        tmp_grad = np.gradient(grad_k)
        for (l, grad_kl) in enumerate(tmp_grad):
            hessian[k, l, :, :] = grad_kl
    return hessian


# Diagonalization functions

def fluxnet_operators(D, N):
    """
    accepts: [D] boundary of the phi grid, [N] number of points of the grid
    returns: operators for Hamiltonian in phaselike basis
    """

    grid_pts = np.linspace(-D, D, N, endpoint=True, dtype='float64')
    dD = 2 * D / (N - 1)

    # Identity

    p0 = ssp.identity(N, format='csc', dtype='float64')

    # Phi

    p1 = ssp.diags(grid_pts, offsets=0, shape=(N, N), format='csc',
                   dtype='float64')

    # Phi**2

    p2 = ssp.diags(grid_pts ** 2, offsets=0, shape=(N, N), format='csc'
                   , dtype='float64')

    # Sin(Phi)

    si = ssp.diags(np.sin(grid_pts), offsets=0, shape=(N, N),
                   format='csc', dtype='float64')

    # Cos(Phi)

    co = ssp.diags(np.cos(grid_pts), offsets=0, shape=(N, N),
                   format='csc', dtype='float64')

    # First derivative

    d1 = 1 / 2 / dD * ssp.diags([-1, 1], offsets=[-1, 1], shape=(N, N),
                                format='csc', dtype='float64')

    # Second derivative

    d2 = 1 / dD ** 2 * ssp.diags([1, -2, 1], offsets=[-1, 0, 1],
                                 shape=(N, N), format='csc',
                                 dtype='float64')
    return [
        p0,
        p1,
        p2,
        si,
        co,
        d1,
        d2,
        ]


def H_inductive(
    phi_ext,
    EC,
    EJ,
    EL,
    ):
    """
    accepts: [phi_ext] phase associated to the external flux, 
    [r] is the ratio of loop sizes, [EC, EJ, EL] parameters
    returns: Hamiltonian in phaselike basis
    non fatso version
    oops. Was this supposed to be the whole energy?
    """

    (
        p0,
        p1,
        p2,
        si,
        co,
        d1,
        d2,
        ) = fluxnet_operators(max_grid, N_grid)

    return -4 * EC * d2 + +EL / 2 * p2 - EJ * co * np.cos(phi_ext) - EJ \
        * si * np.sin(phi_ext)


def H_potential(phi_ext, EJ, EL):
    """
    accepts: [phi_ext] phase associated to the external flux, [r] is the 
    ratio of loop sizes, [EC, EJ, EL] parameters
    returns: Hamiltonian in phaselike basis
    
    Non fatso version
    """

    grid_pts = np.linspace(-max_grid, max_grid, N_grid, endpoint=True,
                           dtype='float64')

    return -EJ * np.cos(grid_pts - phi_ext) + +EL / 2 * grid_pts ** 2


def diagonalize_new(operator, number_of_ekets):
    """
    returns [evals, ekets]: sorted in energy 
    """

    (evals, ekets) = scipy.sparse.linalg.eigsh(operator,
            k=number_of_ekets, which='SA', tol=tol)
    indx = evals.argsort()
    evals_s = np.sort(evals)
    evals_s = evals_s
    ekets_s = np.zeros(ekets.shape, dtype='float64')
    for i in range(number_of_ekets):
        ekets_s[:, i] = ekets[:, indx[i]]
    return (evals_s, ekets_s)


def flux_sweep(
    phi_ext_s,
    EC,
    EJ,
    EL,
    ):
    """
    returns [ evals_s, ekets_s]: sorted in energy 
    """

    evals_s = []
    ekets_s = []
    for phi_ext in phi_ext_s:
        H_op = H_inductive(phi_ext, EC, EJ, EL)
        (temp_evals, temp_ekets) = diagonalize_new(H_op, keig)
        evals_s.append(temp_evals)
        ekets_s.append(temp_ekets)
    evals_s = np.array(evals_s)
    ekets_s = np.array(ekets_s)
    return (evals_s, ekets_s)


def g_dipole_calc(e_vals, e_kets):
    """
    e_vals :  keigs dim array containing eigenvalues. ( don't really need it)
    e_kets :  Nxkeigs dimensional array containing eigenvectors in each column 
    """

    keigs = len(e_vals)
    g_matrix = np.zeros((keigs, keigs))
    delta_matrix = np.zeros((keigs, keigs))
    chi_matrix = np.zeros((keigs, keigs))
    chi_s = np.zeros(keigs)

    D = max_grid
    N = N_grid
    grid_pts = np.linspace(-D, D, N, endpoint=True, dtype='float64')
    dD = grid_pts[1] - grid_pts[0]

    # First derivative

    d1 = 1 / 2 / dD * ssp.diags([-1, 1], offsets=[-1, 1], shape=(N, N),
                                format='csc', dtype='float64')
    for i in np.arange(keigs):
        for j in np.arange(keigs):
            g_matrix[i, j] = g_phi * np.vdot(e_kets[:, i],
                    d1.dot(e_kets[:, j]))
            delta_matrix[i, j] = e_vals[i] - e_vals[j] - f_res
            chi_matrix[i, j] = g_matrix[i, j] ** 2 / delta_matrix[i, j]
    chi_s = np.sum(chi_matrix, 1) - np.sum(chi_matrix, 0)
    return (g_matrix, chi_s)


def g_phi_calc(phi_ext, e_vals, e_kets):
    """
    e_vals :  keigs dim array containing eigenvalues. ( don't really need it)
    e_kets :  Nxkeigs dimensional array containing eigenvectors in each column 
    """

    keigs = len(e_vals)
    g_matrix = np.zeros((keigs, keigs))
    D = max_grid
    N = N_grid
    grid_pts = np.linspace(-D, D, N, endpoint=True, dtype='float64')

    # First derivative

    p1 = ssp.diags(grid_pts - phi_ext, offsets=0, shape=(N, N),
                   format='csc', dtype='float64')
    for i in np.arange(keigs):
        for j in np.arange(keigs):
            g_matrix[i, j] = np.vdot(e_kets[:, i], p1.dot(e_kets[:, j]))
    return g_matrix


def g_quasi_calc(phi_ext, e_vals, e_kets):
    """
    e_vals :  keigs dim array containing eigenvalues. ( don't really need it)
    e_kets :  Nxkeigs dimensional array containing eigenvectors in each column 
    """

    keigs = len(e_vals)
    g_matrix = np.zeros((keigs, keigs))
    D = max_grid
    N = N_grid
    grid_pts = np.linspace(-D, D, N, endpoint=True, dtype='float64')

    # First derivative

    si2 = ssp.diags(np.sin((grid_pts - phi_ext) / 2), offsets=0,
                    shape=(N, N), format='csc', dtype='float64')
    for i in np.arange(keigs):
        for j in np.arange(keigs):
            g_matrix[i, j] = np.vdot(e_kets[:, i], si2.dot(e_kets[:,
                    j]))
    return g_matrix


def g_dipole_flux_sweep(  # from Pranav: doesn't say much about T1
    phi_ext_s,
    EC,
    EJ,
    EL,
    ):
    """
    returns [ evals_s]: sorted in energy 
    """

    evals_s = []
    g_flux = []
    for phi_ext in phi_ext_s:
        H_op = H_inductive(phi_ext, EC, EJ, EL)
        (temp_evals, temp_ekets) = diagonalize_new(H_op, keig)
        (temp_g_mat, temp_chi) = g_dipole_calc(temp_evals, temp_ekets)
        evals_s.append(temp_evals)
        g_flux.append([temp_g_mat[0, 1], temp_g_mat[0, 2], temp_g_mat[1, 2]])
        #g_flux.append([temp_g_mat[0, 1])
    evals_s = np.array(evals_s)
    return (evals_s, g_flux)


def g_phi_flux_sweep(
    phi_ext_s,
    EC,
    EJ,
    EL,
    ):
    """
    returns [ evals_s]: sorted in energy 
    """

    evals_s = []
    g_flux = []
    for phi_ext in phi_ext_s:
        H_op = H_inductive(phi_ext, EC, EJ, EL)
        (temp_evals, temp_ekets) = diagonalize_new(H_op, keig)
        temp_g_mat = g_phi_calc(phi_ext, temp_evals, temp_ekets)
        evals_s.append(temp_evals)
        g_flux.append([temp_g_mat[0, 1], temp_g_mat[0, 2], temp_g_mat[1, 2]])
    evals_s = np.array(evals_s)
    return (evals_s, np.array(g_flux))


def g_quasi_flux_sweep(
    phi_ext_s,
    EC,
    EJ,
    EL,
    ):
    """
    returns [ evals_s]: sorted in energy 
    """

    evals_s = []
    g_flux = []
    for phi_ext in phi_ext_s:
        H_op = H_inductive(phi_ext, EC, EJ, EL)
        (temp_evals, temp_ekets) = diagonalize_new(H_op, keig)
        temp_g_mat = g_quasi_calc(phi_ext, temp_evals, temp_ekets)
        evals_s.append(temp_evals)
        g_flux.append(temp_g_mat[0, 1])
    evals_s = np.array(evals_s)
    return (evals_s, g_flux)


def double_flux_sweep(
    phi1_ext_s,
    phi2_ext_s,
    EC,
    EJ,
    EL,
    ):
    """
    returns [ evals_s]: sorted in energy 
    """

    N_phi1 = phi1_ext_s.shape[0]
    N_phi2 = phi2_ext_s.shape[0]
    evals_s = np.zeros((N_phi1, N_phi2, keig))
    for i in range(N_phi1):
        for j in range(N_phi2):
            phi_ext1 = phi1_ext_s[i]
            phi_ext2 = phi2_ext_s[j]
            r = phi_ext1 / phi_ext2
            H_op = H_inductive(phi_ext2, EC, EJ, EL)
            (temp_evals, temp_ekets) = diagonalize_new(H_op, keig)
            evals_s[i, j, :] = temp_evals
    return evals_s

###################################################################
###### making Raman functions
def RamanVal(
    EC,
    EJ,
    EL,
    cavfreq,
    phi_ext,
    numstart=0,
    drive_freq = 0.50,
    delta_lil = 0.05,
    ):
    keig = 10  # number of levels to be summed over

    # ## finding the raman transition

    H_op = H_inductive(phi_ext, EC, EJ, EL)
    (evals_s, ekets_s) = diagonalize_new(H_op, keig)
    evals_s = evals_s - evals_s[numstart]

    # finding the dipole matrix

    (g_matrix, chi_s) = g_dipole_calc(evals_s, ekets_s)

    deltas = evals_s - cavfreq

    gconst = cavfreq * beta_phi * (2 / np.pi)

    # ### making drive consistant for all phi at 50 MHz

    drive = drive_freq / g_matrix[numstart, 2]

    # ###

    gvals = gconst * g_matrix[:, 1]
    omgvals = g_matrix[0, :] * drive

    ramanvals = np.zeros(len(evals_s))
    for i in range(1, len(evals_s)):
        ramanvals[i] = gvals[i] * omgvals[i] /((deltas[i] - delta_lil)*2)

    OmgRaman = np.sum(ramanvals) 

    return OmgRaman


##########################################################################
def ChiShift(
    EC,
    EJ,
    EL,
    cavfreq,
    phi_ext,
    level = 0,
    ):

    keig = 10  # number of levels to be summed over
    ChiNaN = False
    # ## finding the energies  
    H_op = H_inductive(phi_ext, EC, EJ, EL)
    (evals_s, ekets_s) = diagonalize_new(H_op, keig)
    evals_s = evals_s - evals_s[0]

    # finding the dipole matrix
    (g_matrix, chi_s) = g_dipole_calc(evals_s, ekets_s)
    gconst = cavfreq * beta_phi * (2 / np.pi)
    gvals = gconst * g_matrix[:, level]
    Chivals = np.zeros(len(evals_s))

    for i in range(0, len(evals_s)):
        Chi1 = np.abs(gvals[i])**2*(1/(evals_s[level] - evals_s[i] - cavfreq))
        Chi2 = np.abs(gvals[i])**2*(1/(evals_s[i] - evals_s[level] - cavfreq))
        Chivals[i] = Chi1 - Chi2 
        ### get right of self case
        if i == level:
            Chivals[i] = 0
        #### looking for non dispersive components
        delt1 = (1/(evals_s[level] - evals_s[i] - cavfreq))
        delt2 = (1/(evals_s[i] - evals_s[level] - cavfreq))
        frac = 1/(.33)
        
        if np.abs(delt1) > frac or np.abs(delt2) > frac:
            ChiNaN = True
            break

    ChiShift = (np.sum(Chivals))

    if ChiNaN:
        ChiShift = float('NaN')

    return ChiShift

###########################################################################

def gfourthVal(
    EC,
    EJ,
    EL,
    cavfreq,
    phi_ext,
    numstart=0,
    ):
    keig2 = 10  # number of levels to be summed over

    # ### finding the g4 value

    g4List = np.ones(keig2 ** 3)
    H_op = H_inductive(phi_ext, EC, EJ, EL)
    (evals_s, ekets_s) = diagonalize_new(H_op, keig2)
    evals_s = evals_s - evals_s[0]  # setting zero energy

    # #### finding the dipole matrix

    (g_matrix, chi_s) = g_dipole_calc(evals_s, ekets_s)
    gconst = cavfreq * beta_phi * (2 / np.pi)
    g_matrix = g_matrix * gconst

    index = 0
    numphot = 0
    for ind1 in range(0, keig2):

        gind01 = g_matrix[numstart, ind1]
        E1 = evals_s[ind1] - cavfreq

        for ind2 in range(0, keig2):

            gind12 = g_matrix[ind1, ind2]
            numphot = -1 - np.sign(ind2 - ind1)
            E2 = evals_s[ind2] + numphot * cavfreq

            if ind2 == 0:
                E2 = E1

            for ind3 in range(0, keig2):

                gind23 = g_matrix[ind2, ind3]
                gind30 = g_matrix[ind3, numstart]
                numphot = -1 - np.sign(ind2 - ind1) - np.sign(ind3
                        - ind2)
                E3 = evals_s[ind3] + numphot * cavfreq
                g4List[index] = gind01 * gind12 * gind23 * gind30 / (E1
                        * E2 * E3)

                numphot = numphot + 1
                if ind1 == ind2 or ind2 == ind3 or ind1 == numstart \
                    or ind3 == numstart or numphot != 0:
                    g4List[index] = 0

                if ind1 == ind3 and ind2 == 0:
                    g4List[index] = 0

                if ind1 == 1 or ind2 == 1 or ind3 == 1:
                    pass
                else:
                    g4List[index] = 0

                index += 1

    # #######

    g01 = g_matrix[0, 1]
    g10 = g_matrix[1, 0]
    g1shift = g10 * g01 / (evals_s[1] - cavfreq)
    g4List = np.append(g4List, g1shift)

    # ######

    gfourth = np.sum(g4List)

    return gfourth

###########################################################################
