'''
This module contains the refract function which implements Snell's Law of refraction at a surface.
'''
import numpy as np
def refract(direc, normal, n_1, n_2):
    '''
    This function returns the refracted ray direction.
    
    Args:
        direc (np.ndarray): the incident direction of the ray
        normal (np.ndarray): the surface noraml
        n_1 (float): refractive index of one side
        n_2 (float): refractive index of other side
    
    Returns:
        np.ndarray: Normalised Numpy array of the refracted ray's direction.
    '''
    #normalise the direction and normal vectors first:
    direc /= np.sqrt(np.dot(direc,direc))
    normal /= np.sqrt(np.dot(normal,normal))

    #incident angle:
    if direc[-1]/normal[-1] < 0:
        cos_theta_1 = -np.dot(direc,normal)
    else:
        cos_theta_1 = np.dot(direc, normal)
        normal*=-1
    sin_theta_1 = np.sqrt(1 - cos_theta_1**2)

    #refracted angle:
    sin_theta_2 = (n_1/n_2) * sin_theta_1
    cos_theta_2 = np.sqrt(1 - sin_theta_2**2)

    #refracted ray direction:
    refracted_direc = ((n_1 / n_2) * np.array(direc)) + ((cos_theta_1 * (n_1 / n_2)) - cos_theta_2) * np.array(normal)
    refracted_direc /= np.sqrt(np.dot(refracted_direc,refracted_direc))

    if sin_theta_1 > (n_2/n_1):
        return None
    else:
        return np.array(refracted_direc)