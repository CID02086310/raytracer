"""
This module contains the Ray class which represents an optical ray.
"""
import matplotlib.pyplot as plt
import numpy as np
class Ray:
    """"
    This class represents an optical ray. 
    Each point as well as the normalised direction is a 3D vector in Cartesian coordinates.
    """
    def __init__(self, pos = np.array([0.,0.,0.]), direc = [0.,0.,1.]):
        '''
        initialises the position, direction and the list of points with the starting position

        Args:
            pos (np.ndarray): initial position
            direc (np.ndarray): initial direction
        
        Raises:
            TypeError: If pos or direc is not a list or a numpy array.
            ValueError: If pos or direc does not have exactly 3 elements.
        '''
        if not isinstance(pos, (list, np.ndarray)):
            raise TypeError("pos must be a list or a numpy array")
        if len(pos) != 3:
            raise ValueError("pos must have exactly 3 elements")
        if not isinstance(direc, (list, np.ndarray)):
            raise TypeError("direc must be a list or a numpy array")
        if len(direc) != 3:
            raise ValueError("direc must have exactly 3 elements")

        self.__pos = np.array(pos)
        self.__direc = np.array(direc)
        self.__direc /= np.sqrt(np.dot(np.array(direc), np.array(direc)))

        self.__points = [self.__pos.copy()]

    def pos(self): #getter for position
        '''
        Returns position
        
        Returns:
            np.ndarray: The current position vector
        '''
        return np.array(self.__pos)
    def direc(self): # getter for direction
        '''
        Returns direction

        Returns: 
            np.ndarray: The current direction vector
        '''
        return np.array(self.__direc)

    def append(self, pos, direc):
        '''
        Appends a new point and direction to the ray

        Args: 
            pos (np.ndarray): the new position
            direc (np.ndarray): the new direction

        Returns:
            np.ndarray: The new position vector and the new direction vector

        Raises:
            TypeError: If pos or direc is not a list or a numpy array.
            ValueError: If pos or direc does not have exactly 3 elements.
        '''
        if not isinstance(pos, (list, np.ndarray)):
            raise TypeError("pos must be a list or a numpy array")
        if len(pos) != 3:
            raise ValueError("pos must have exactly 3 elements")

        if not isinstance(direc, (list, np.ndarray)):
            raise TypeError("direc must be a list or a numpy array")
        if len(direc) != 3:
            raise ValueError("direc must have exactly 3 elements")
        self.__pos = np.array(pos)
        self.__direc = np.array(direc)
        self.__direc /= np.sqrt(np.dot(np.array(direc), np.array(direc)))
        self.__points.append(self.__pos.copy())
        return self
    def vertices(self):
        '''
        Returns all the points along the ray in the form of a list.
        
        Returns:
            np.ndarray: array containing the initial and all appended position vectors   
        '''
        return self.__points
ray1= Ray()
ray1.append([10., 12., 14.], [11., 13., 15.])

class RayBundle:
    '''
    Represents a circular (transverse to the optical axis) bundle of Ray objects.'''
    def __init__(self, rmax = 5., nrings = 5, multi = 6):
        '''
        initialises the values of rmax, nrings and multi in the RayBundle class.
        
        Args:
            rmax (float): the maximum radius the circular bundle of ray objects extends to
            nrings (float): the number of concentric rings to produce in the circular bundle
            multi (float): a multiplier which indicates how many points should be '''
        self.__rmax = rmax
        self.__nrings = nrings
        self.__multi = multi
        self.__rays  =self.rtrings()

    def rtrings(self):
        '''
        Constructs the bundle of ray objects
        
        Returns:
            list: a list containing all the rays in the bundle'''
        rays = []
        rays.append(Ray(pos = np.array([0., 0., 0.]), direc =[0., 0., 1.] ))

        # Calculate the step size for radius
        radius_step = self.__rmax / self.__nrings
        for i in range(self.__nrings):
            radius = radius_step * (i + 1)
            num_points = (i + 1) * self.__multi
            for j in range(num_points):
                theta = 2 * np.pi * j / num_points
                x = radius * np.cos(theta)
                y = radius * np.sin(theta)
                rays.append(Ray(pos = np.array([x, y, 0.]), direc = [0., 0., 1.]))
        return rays
    def propagate_bundle(self, elements):
        '''
        Takes a list of Optical Element objects and sequentially propagates each ray through them.

        Args: 
            elements (list): a list of optical element objects
        '''
        for ray in self.__rays:
            for element in elements:
                element.propagate_ray(ray)
    def track_plot(self):
        '''
        Creates a plot of the path of each ray in the bundle, returns the matplotlib figure object.

        Returns: 
            matplotlib.pyplot: a figure, displaying the path of each ray in the bundle.
        '''
        fig = plt.figure()
        ax = fig.add_subplot(projection = "3d") #this is to plot the graph in 3d
        for ray in self.__rays:
            vertices =ray.vertices()
            vertices = np.array(vertices)
            x = []
            y = []
            z = []
            for intercept in vertices:
                x.append(intercept[0])
                y.append(intercept[1])
                z.append(intercept[2])
            plt.plot(x,y,z)
        plt.xlabel("z position")
        plt.ylabel("x position")
        plt.title("Plot of each ray in the bundle")

        return fig

    def rms(self):
        '''
        
        This method calculates and returns the RMS size of the bundle.
        
        Returns:
            float: the RMS size of the bundle.
        '''
        x = [ray.pos()[0] for ray in self.__rays]
        y = [ray.pos()[1] for ray in self.__rays]

        x = np.array(x)
        y = np.array(y)

        rms= np.sqrt(np.mean(x**2  + y**2))
        return  rms

    def spot_plot(self):
        '''
        Plots the x-y position of each ray in the bundle at the rays' current positions
        
        Returns:
            matplotlib.pyplot: the x-y pos of each ray in the bundle at the rays' current positions.
        '''
        fig = plt.figure()
        x = [ray.pos()[0] for ray in self.__rays]
        y = [ray.pos()[1] for ray in self.__rays]
        plt.scatter(x,y)
        plt.xlabel("x position")
        plt.ylabel("y position")
        plt.title("The x-y position of each ray in the bundle")
        return fig