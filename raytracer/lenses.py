'''The lenses module contains the Convex-plano and Plano-convex lenses'''
import numpy as np
from raytracer.elements import OpticalElement
from raytracer.elements import SphericalRefraction
class PlanoConvex(OpticalElement):
    '''
    This class represents Plano Convex lenses'''
    def __init__(self, z_0, curvature1, curvature2, n_inside, n_outside, thickness, aperture):
        '''
        Initialises the intercept, the curvatures of the surface, n_inside, n_outside and the thickness of the lens
        
        Args:
            z_0 (float): the intercept of the left most surface with the z-axis
            curvature1 (float): the curvature of the left surface
            curvature2 (float): the curvature of the right surface
            n_inside (float): the refractive index of the inside of the lens
            n_outside (float): the refractive index of the outside of the lens
            thickness (float): the thickness of the lens on the optical (z) axis
            aperture (float): the maximum extent of the surface from the optical axis
        '''
        self.__z_0 = z_0
        self.__curvature1 = curvature1
        self.__curvature2 = curvature2
        self.__n_inside = n_inside
        self.__n_outside = n_outside
        self.__thickness = thickness
        self.__aperture = aperture

        if self.__curvature1 == 0:
            self.__curvature1 = 1e-12

        if self.__curvature2 == 0:
            self.__curvature2 = -1e-12

        self.__side1 = SphericalRefraction(self.__z_0, self.__aperture, self.__curvature1,
                                           self.__n_outside, self.__n_inside)
        self.__side2 = SphericalRefraction(self.__z_0 + self.__thickness, self.__aperture,
                                           self.__curvature2, self.__n_inside, self.__n_outside)

    def side1(self):
        '''
        Returns the object for the first side of the lens

        Returns:
            SphericalRefraction: the first side of the lens.
        '''
        return self.__side1
    def side2(self):
        '''
        Returns the obeject for the second side of the lens

        Returns:
            SphericalRefraction: the second side of the lens.
        '''
        return self.__side2

    def z_0(self):
        '''
        Returns the value of the intercept of the left most surface with the z-axis
        
        Returns:
            float: the value of the intercpet, z_0
        '''
        return self.__z_0
    def curvature1(self):
        '''
        Returns the value of the curvature of the left surface
        
        Returns:
            float: curvature of the left surface
        '''
        return self.__curvature1

    def curvature2(self):
        '''
        Returns the value of the curvature of the right surface
        
        Returns:
            float: curvature of the right surface
        '''
        return self.__curvature2
    def n_inside(self):
        '''
        Returns the value of the refractive index of the inside of the lens
        
        Returns:
            float: the refractive index of the inside of the lens
        '''
        return self.__n_inside
    def n_outside(self):
        '''
        Returns the value of the refractive index of the outside of the lens.
        
        Returns: 
            float: the refractive index of the outside of the lens
        '''
        return self.__n_outside
    def thickness(self):
        '''
        Returns the value of the thickness of the lens on the optical axis
        
        Returns:
            float: the thickness of the lens
        '''
        return self.__thickness
    def aperture(self):
        '''
        Returns the value of the aperture of the lens
        
        Returns: 
            float: the aperture
        '''
        return self.__aperture
    def focal_point(self):
        '''
        This method calculates and returns the focal point of the lens.
         
        Returns:
            float: the focal point of the lens
        '''

        if self.__curvature1 == 1e-12:
            R = 1/self.__curvature2
            focal_length = (R) / (self.__n_outside - self.__n_inside)
            focal_point_pc = self.__z_0 + focal_length + self.__thickness
            return focal_point_pc
        if self.__curvature2 == -1e-12:
            R = 1/self.__curvature1
            focal_point_cp = self.__z_0 + self.__thickness + (
                (R / (self.__n_inside - self.__n_outside)) - (self.__thickness / self.__n_inside))
            return focal_point_cp
        
class ConvexPlano(OpticalElement):
    """
    This class represents Convex Plano lenses.
    """

    def __init__(self, z_0, curvature, n_inside, n_outside, thickness, aperture):
        """
        Initializes the intercept, curvature, n_inside, n_outside, & the thickness.

        Args:
            z_0 (float): The intercept of the left-most surface with the z-axis.
            curvature (float): The curvature of the convex surface.
            n_inside (float): The refractive index of the inside of the lens.
            n_outside (float): The refractive index of the outside of the lens.
            thickness (float): The thickness of the lens on the optical (z) axis.
            aperture (float): The maximum extent of the surface from the optical axis.
        """
        self.__z_0 = z_0
        self.__curvature = curvature
        self.__n_inside = n_inside
        self.__n_outside = n_outside
        self.__thickness = thickness
        self.__aperture = aperture

        if self.__curvature == 0:
            self.__curvature = 1e-12

        self.__side1 = SphericalRefraction(self.__z_0, self.__aperture, self.__curvature, 
                                           self.__n_outside, self.__n_inside)
        self.__side2 = SphericalRefraction(self.__z_0 + self.__thickness, self.__aperture,
                                           -1e-12, self.__n_inside, self.__n_outside)

    def side1(self):
        """
        Returns the object for the first side of the lens.

        Returns:
            SphericalRefraction: The first side of the lens.
        """
        return self.__side1

    def side2(self):
        """
        Returns the object for the second side of the lens.

        Returns:
            SphericalRefraction: The second side of the lens.
        """
        return self.__side2

    def z_0(self):
        """
        Returns the value of the intercept of the left-most surface with the z-axis.

        Returns:
            float: The value of the intercept, z_0.
        """
        return self.__z_0

    def curvature(self):
        """
        Returns the value of the curvature of the convex surface.

        Returns:
            float: The curvature of the convex surface.
        """
        return self.__curvature

    def n_inside(self):
        """
        Returns the value of the refractive index of the inside of the lens.

        Returns:
            float: The refractive index of the inside of the lens.
        """
        return self.__n_inside

    def n_outside(self):
        """
        Returns the value of the refractive index of the outside of the lens.

        Returns:
            float: The refractive index of the outside of the lens.
        """
        return self.__n_outside

    def thickness(self):
        """
        Returns the value of the thickness of the lens on the optical axis.

        Returns:
            float: The thickness of the lens.
        """
        return self.__thickness

    def aperture(self):
        """
        Returns the value of the aperture of the lens.

        Returns:
            float: The aperture.
        """
        return self.__aperture

    def focal_point(self):
        """
        This method calculates and returns the focal point of the lens.

        Returns:
            float: The focal point of the lens.
        """
        R = 1 / self.__curvature
        focal_length = R / (self.__n_outside - self.__n_inside)
        focal_point = self.__z_0 + focal_length + self.__thickness
        return focal_point