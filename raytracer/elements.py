'''
This module contains the OpticalElement class and SphericalRefraction class
Classes:
    OpticalElement: Base class for optical elements.
    SphericalRefraction: Represents a spherical refracting surface.
    OutputPlane: Represents the plane where rays are propagated to but do not refract.
'''
import numpy as np
from raytracer.physics import refract

class OpticalElement:
    '''
    The OpticalElement class is used to represent an optical system.
    Methods:
        intercept(ray): Calculate where a ray intercepts the element.
        propagate_ray(ray): Propagate a ray through the element.
    '''

    def intercept(self,ray):
        '''
        Takes a ray object as an argument & calculates where its path intercepts with the element.

        Args:
            ray (Ray): an object, representing the incoming ray

        Raises:
            NotImplementedError: If the method is not implemented in a derived class.
        '''
        raise NotImplementedError('intercept() needs to be implemented in derived classes')

    def propagate_ray(self, ray):
        """
        This method takes a ray object as an argument and will propagate it through the element.

        Args:
            ray (Ray): An object representing the incoming ray.

        Raises:
            NotImplementedError: If the method is not implemented in a derived class.
        """
        raise NotImplementedError('propagate_ray() needs to be implemented in derived classes')

class Sphere(OpticalElement):
    '''
    Represents a Sphere.
    '''
    def __init__(self, z_0, aperture, curvature):
        '''
        initialises the SphericalRefraction witht he given parameters.
        
        Args:
            z_0 (float): the intercept of the surface with the z-axis.
            aperture (float): the maximum extent of the surface from the optical axis.
            curvature (float):  the curvature of the surface (magnitude 1/(radius of curvature)).
        '''
        self.__z_0 = z_0
        self.__aperture = aperture
        self.__curvature = curvature
    def z_0(self):
        '''
        Returns the intercept of the surface with the z-axis.
        
        Returns: 
            float: the inercept of the surface with the z-axis
        '''
        return self.__z_0

    def aperture(self):
        '''
        Returns the aperture (maximum extent of the surface from the optical axis)
        
        Returns: 
            float: the aperture of the surface.
        '''
        return self.__aperture
    def curvature(self):
        '''
        Returns the curvature of the surface.
        
        Returns: 
            float: the curvature
        '''
        return self.__curvature

    def intercept(self, ray):
        '''
        Finds Q, the point of intersection of the incoming ray from point P, the source of the ray.
         
        Args: 
            ray (Ray): an object, representing the incoming ray
             
        Returns: 
            None: If no valid intercept
            np.ndarray: numpy array of the intersection point
        '''
        P = ray.pos()
        k = ray.direc()
        R = 1/self.curvature() # radius of curvature
        O = np.array([0,0,self.z_0() + R])

        r = P - O
        r_dot_k = np.dot(r,k) 
        r_squared = np.dot(r,r) 
        discriminant = (r_dot_k * r_dot_k) - (r_squared - R**2)
        l1 = -(r_dot_k) + np.sqrt(discriminant)
        l2 = -(r_dot_k) - np.sqrt(discriminant)

        if discriminant < 0:
            return None
        if l1 <= 0 and l2 <= 0:
            return None
        if l2 < 0 and l1 > 0:
            l = l1
        if l1 < 0 and l2 > 0:
            l = l2
        if l1 > 0 and l2 > 0:
            l = min(l1,l2)
        Q = P + (l * k)
        if Q[0] >self.aperture() or Q[1] > self.aperture():
            return None
        return Q

class SphericalRefraction(Sphere):
    '''
    Represents a spherical refracting surface.'''

    def __init__(self, z_0, aperture, curvature, n_1, n_2):
        '''
        initialises the SphericalRefraction witht he given parameters.
        
        Args:
            z_0 (float): the intercept of the surface with the z-axis.
            aperture (float): the maximum extent of the surface from the optical axis.
            curvature (float):  the curvature of the surface (magnitude 1/(radius of curvature)).
            n_1 (float): the refractive index on the left of the mirror.
            n_2 (float): the refractive index on the right of the mirror.
        '''
        self.__z_0 = z_0
        self.__aperture = aperture
        self.__curvature = curvature
        self.__n_1 = n_1
        self.__n_2 = n_2

    def z_0(self):
        '''
        Returns the intercept of the surface with the z-axis.
        
        Returns: 
            float: the inercept of the surface with the z-axis
        '''
        return self.__z_0

    def aperture(self):
        '''
        Returns the aperture (maximum extent of the surface from the optical axis)
        
        Returns: 
            float: the aperture of the surface.
        '''
        return self.__aperture
    def curvature(self):
        '''
        Returns the curvature of the surface.
        
        Returns: 
            float: the curvature
        '''
        return self.__curvature
    def n_1(self):
        '''
        Returns the refractive index on the left of the surface.

        Returns: 
            float: left side refractive index
        '''
        return self.__n_1
    def n_2(self):
        '''
        Returns the refractive index on the right of the surface.
        
        Returns: 
            float: right side refractive index
        '''
        return self.__n_2
    def propagate_ray(self, ray):
        '''
        propagates the ray object through the surface of the optical element
        
        Args:
            ray (Ray): an object, representing an incoming ray
            
        Returns:
            np.ndarray: the new position, new direction
        '''
        R = 1/self.__curvature # radius of curvature
        O = np.array([0,0,self.__z_0 + R]) 
        new_pos =self.intercept(ray)
        if new_pos is None:
            print("No valid intercept is found.")
            return None
        normal = np.array(new_pos) - np.array(O)

        new_direc = refract(ray.direc(), normal, self.__n_1, self.__n_2)
        if new_direc is None:
            print("Total Internal Reflaction occured.")
            return None
        ray.append(pos = new_pos, direc = new_direc)

    def focal_point(self):
        '''
        This method calculates and returns the focal point of the surface.
        
        Returns:
            float: the focal point of the surface'''
        R = 1/self.__curvature # radius of curvature
        focal_length = (self.__n_2 * R) / (self.__n_2 - self.__n_1)
        focal_point = self.__z_0 + focal_length

        return focal_point
class Plane(OpticalElement):
    '''
    Represents a plane.
    '''
    def __init__(self,z_0):
        '''
        Initialises the output plane and the intercept of the surface with the z-axis
        
        Args:
            z_0 (float): the intercept of the surface with the z-axis'''
        self.__z_0 = z_0
    def z_0(self):
        '''
        Returns the intercept of the surface with the z-axis.
        
        Returns: 
            float: the inercept of the surface with the z-axis
        '''
        return self.__z_0
    def intercept(self,ray):
        '''
        Find the intersection of the ray with the output plane.
        
        Args:
            ray (Ray): an object, representing an incoming ray
            
        Returns:
            np.ndarray: intercept, or None if the ray is parallel to the output plane.
        '''
        P = ray.pos()
        k = ray.direc()

        if k[2] == 0:
            return None # The ray is parallel to the output plane
        l = (self.z_0() - P[2]) / k[2] 
        intercept = P + l*k

        return intercept

class OutputPlane(Plane):
    '''
    Rays are propagated to the output plane but their direction remains unchanged'''
    def __init__(self,z_0):
        '''
        Initialises the output plane and the intercept of the surface with the z-axis
        
        Args:
            z_0 (float): the intercept of the surface with the z-axis'''
        self.__z_0 = z_0

    def z_0(self):
        '''
        Returns the intercept of the surface with the z-axis.
        
        Returns: 
            float: the inercept of the surface with the z-axis
        '''
        return self.__z_0

    def propagate_ray(self, ray):
        '''
        Propagates the ray to the output plane and appends the new position of the ray.
        
        Args:
            ray (Ray): an object, representing an incoming ray

        Returns:
            aNone: desciption
        '''
        new_pos = self.intercept(ray)

        if new_pos is None:
            print("No valid intercept is found.")
            return None
        new_direc = ray.direc() #The direction of the ray remains unchanged

        ray.append(new_pos, new_direc)
