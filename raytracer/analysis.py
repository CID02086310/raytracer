"""Analysis module."""
import matplotlib.pyplot as plt
import numpy as np
from raytracer.elements import SphericalRefraction
from raytracer.elements import OutputPlane
from raytracer.rays import Ray
from raytracer.rays import RayBundle
from raytracer.lenses import PlanoConvex


def task8():
    """
    Task 8.

    In this function you should check your propagate_ray function properly
    finds the correct intercept and correctly refracts a ray. Don't forget
    to check that the correct values are appended to your Ray object.
    """
    ray1 = Ray(pos = np.array([0.,0.,0.]), direc= [0.,0.,1.]) #original ray: along the axis
    ray2 = Ray(pos = np.array([0.,2.,0.]), direc= [0.,0.,1.])
    ray3 = Ray(pos = np.array([1.,1.,0.]), direc= [0.,0.,1.]) #off-axis ray

    rays= [ray1,ray2,ray3]
    test_pos = []
    test_direc = []
    sr = SphericalRefraction(z_0 = 10, aperture = 5, curvature = 0.2, n_1 = 1, n_2 = 1.5)

    for ray in rays:
        sr.propagate_ray(ray) # pass the ray through the optical element (sr)
        test_pos.append(ray.pos())
        test_direc.append(ray.direc())
        print(ray.pos())
        print(ray.direc())
        # assert np.allclose(np.array(test_pos),np.array(expected_pos))
        # assert np.allclose(test_direc, expected_direc)
    print(test_pos)
def task10():
    """
    Task 10.

    In this function you should create Ray objects with the given initial positions.
    These rays should be propagated through the surface, up to the output plane.
    You should then plot the tracks of these rays.
    This function should return the matplotlib figure of the ray paths.

    Returns:
        Figure: the ray path plot.
    """
    sr = SphericalRefraction(z_0 = 100, aperture = 5, curvature = 0.03, n_1 = 1, n_2 = 1.5)
    op = OutputPlane(250)
    ray1= Ray(pos = np.array([0., 4., 0.]), direc = [0., 0., 1.])
    ray2= Ray(pos = np.array([0., 1., 0.]), direc  =[0., 0., 1.])
    ray3= Ray(pos = np.array([0., 0.2, 0.]), direc = [0., 0., 1.])
    ray4= Ray(pos = np.array([0., 0., 0.]), direc = [0., 0., 1.])
    ray5= Ray(pos = np.array([0., -0.2, 0]), direc = [0., 0., 1.])
    ray6= Ray(pos = np.array([0., -1., 0.]), direc = [0., 0., 1.])
    ray7= Ray(pos = np.array([0., -4., 0.]), direc = [0., 0., 1.])

    initial_rays = [ray1, ray2, ray3, ray4, ray5, ray6, ray7]

    fig = plt.figure()

    for ray in initial_rays:
        sr.propagate_ray(ray) #pass the ray through the refracting sphere first
        op.propagate_ray(ray) #now pass the ray through the output plane
        vertices = ray.vertices()
        vertices = np.array(vertices)
        y = []
        z= []
        for intercept in vertices:
            y.append(intercept[1])
            z.append(intercept[2])
        plt.plot(z,y)
    plt.xlabel("z position")
    plt.ylabel("x position")
    plt.title("Task 10 Plot")
    return  fig
task10()

def task11():
    """
    Task 11.

    In this function you should propagate the three given paraxial rays through the system
    to the output plane and the tracks of these rays should then be plotted.
    This function should return the following items as a tuple in the following order:
    1. the matplotlib figure object for ray paths
    2. the calculated focal point.

    Returns:
        tuple[Figure, float]: the ray path plot and the focal point
    """

    sr = SphericalRefraction(z_0 = 100, aperture = 5, curvature = 0.03, n_1 = 1, n_2 = 1.5)
    op = OutputPlane(250)

    ray1= Ray(pos = np.array([0.1, 0.1, 0.]), direc = [0., 0., 1.])
    ray2= Ray(pos = np.array([0., 0., 0.]), direc  =[0., 0., 1.])
    ray3= Ray(pos = np.array([-0.1, -0.1, 0.]), direc = [0., 0., 1.])
    paraxial_rays = [ray1,ray2,ray3]

    fig = plt.figure()

    for ray in paraxial_rays:
        f = sr.focal_point()
        sr.propagate_ray(ray) #pass the ray through the refracting sphere first
        op.propagate_ray(ray) #now pass the ray through the output plane
        vertices = ray.vertices()
        vertices = np.array(vertices)
        y = []
        z= []
        for intercept in vertices:
            y.append(intercept[1])
            z.append(intercept[2])
        plt.plot(z,y)
    plt.xlabel("z position")
    plt.ylabel("x position")
    plt.title("Task 11 Plot")
    return  fig, f

task11()

def task12():
    """
    Task 12.

    In this function you should create a RayBunble and propagate it to the output plane
    before plotting the tracks of the rays.
    This function should return the matplotlib figure of the track plot.

    Returns:
        Figure: the track plot.
    """
    rb = RayBundle(rmax=5., nrings=5, multi=6)
    sr = SphericalRefraction(z_0 = 100, aperture = 5, curvature = 0.03, n_1 = 1, n_2 = 1.5)
    f = sr.focal_point()
    op = OutputPlane(f)
    elements = [sr,op]

    rb.propagate_bundle(elements) #propagate the ray bundle through the optical elements
    fig = rb.track_plot() #plot the path of each ray in the bundle

    return fig
task12()


def task13():
    """
    Task 13.

    In this function you should again create and propagate a RayBundle to the output plane
    before plotting the spot plot.
    This function should return the following items as a tuple in the following order:
    1. the matplotlib figure object for the spot plot
    2. the simulation RMS

    Returns:
        tuple[Figure, float]: the spot plot and rms
    """
    rb = RayBundle(rmax = 5, nrings = 5, multi =6)
    sr = SphericalRefraction(z_0 = 100, aperture = 5, curvature = 0.03, n_1 = 1, n_2 = 1.5)
    f = sr.focal_point()
    op = OutputPlane(f)
    elements = [sr,op]

    rb.propagate_bundle(elements)
    fig = rb.spot_plot()
    RMS = rb.rms()


    return fig, RMS
task13()


def task14():
    """
    Task 14.

    In this function you will trace a number of RayBundles through the optical system and
    plot the RMS and diffraction scale dependence on input beam radii.
    This function should return the following items as a tuple in the following order:
    1. the matplotlib figure object for the diffraction scale plot
    2. the simulation RMS for input beam radius 2.5
    3. the diffraction scale for input beam radius 2.5

    Returns:
        tuple[Figure, float, float]: the plot, the simulation RMS value, the diffraction scale.
    """
    lamda = 588e-6
    dx = []
    rms  =[]
    fig = plt.figure()
    r_values = np.linspace(0.1, 10, 100)
    for r in r_values:
        rb = RayBundle(rmax = float(r), nrings = 5, multi =6)
        sr = SphericalRefraction(z_0 = 100, aperture = 10, curvature = 0.03, n_1 = 1, n_2 = 1.5)
        focal_length = sr.focal_point() - 100
        f = sr.focal_point()

        op = OutputPlane(f)
        elements = [sr,op]
        rb.propagate_bundle(elements)
        dx.append((lamda * focal_length) /(2*r))
        rms.append(rb.rms())

    plt.plot(r_values, rms, 'r', label = "RMS")
    plt.plot(r_values, dx, 'g', label = "Diffraction Scale")
    plt.xlabel('Input Beam Radius (D)')
    plt.title('RMS and Diffraction Scale vs Input Beam Radius')
    plt.legend()

    rb_25 = RayBundle(rmax = 2.5, nrings = 5, multi =6)
    sr_25 = SphericalRefraction(z_0 = 100, aperture = 5, curvature = 0.03, n_1 = 1, n_2 = 1.5)
    focal_length_25 = sr_25.focal_point() - 100
    f_25 = sr.focal_point()
    op_25 = OutputPlane(f_25)
    elements_25 = [sr_25,op_25]
    rb_25.propagate_bundle(elements_25)
    rms_25 = rb_25.rms()
    dx_25 = (lamda * focal_length_25) / 5

    return fig, rms_25, dx_25

task14()


def task15():
    """
    Task 15.

    This function should return the following items as a tuple in the following order:
    1. the matplotlib figure object for the spot plot for the plano-convex system
    2. the focal point for the plano-convex lens
    3. the matplotlib figure object for the spot plot for the convex-plano system
    4  the focal point for the convex-plano lens


    Returns:
        tuple[Figure, float, Figure, float]: the spotplots & rms for plano-convex and convex-plano.
    """
    z_0 = 100
    curvature = 0.02
    n_glass = 1.5168
    n_air = 1.0
    thickness = 5
    aperture = 50
    pc = PlanoConvex(z_0, 0, -1*curvature, n_glass, n_air, thickness, aperture)
    f_pc = pc.focal_point()
    op_pc = OutputPlane(f_pc)

    rb_pc = RayBundle(rmax = 5., nrings = 5, multi = 6)

    elements_pc = [pc.side1(), pc.side2(), op_pc]
    rb_pc.propagate_bundle(elements_pc)

    fig_pc = rb_pc.spot_plot()
    plt.title("plano-convex spot plot")
    focal_point_pc = pc.focal_point()

    cp = PlanoConvex(z_0, curvature, 0, n_glass, n_air, thickness, aperture)
    f_cp = cp.focal_point()
    op_cp = OutputPlane(f_cp)

    rb_cp = RayBundle(rmax = 5., nrings = 5, multi = 6)
    elements_cp = [cp.side1(), cp.side2(), op_cp]
    rb_cp.propagate_bundle(elements_cp)

    fig_cp = rb_cp.spot_plot()
    plt.title("convex-plano spot plot")
    focal_point_cp = cp.focal_point()

    return fig_pc, focal_point_pc, fig_cp, focal_point_cp

task15()

def task16():
    """
    Task 16.

    This function should return the following items as a tuple in the following order:
    1. the matplotlib figure object for the diffraction scale plot
    2. the RMS for input beam radius 3.5 for the plano-convex system
    3. the RMS for input beam radius 3.5 for the convex-plano system
    4  the diffraction scale for input beam radius 3.5

    Returns:
        tuple[Figure, float, float, float]: 
        the plot, RMS for plano-convex, RMS for convex-plano, diffraction scale.
    """
    lamda = 588e-6
    z_0 = 100
    curvature = 0.02
    n_glass = 1.5168
    n_air = 1.0
    thickness = 5
    aperture = 50

    dx_pc = []
    rms_pc = []

    dx_cp = []
    rms_cp = []

    fig = plt.figure()
    r_values = np.linspace(0.1, 10, 100)
    for r in r_values:
        rb_pc = RayBundle(rmax=float(r), nrings=5, multi=6)
        pc = PlanoConvex(z_0, 0, -1 * curvature, n_glass, n_air, thickness, aperture)
        focal_length_pc = pc.focal_point() - z_0
        f_pc = pc.focal_point()
        op_pc = OutputPlane(f_pc)
        elements_pc = [pc.side1(), pc.side2(), op_pc]
        rb_pc.propagate_bundle(elements_pc)
        dx_pc.append((lamda * focal_length_pc) / (2 * r))
        rms_pc.append(rb_pc.rms())

    for r in r_values:
        rb_cp = RayBundle(rmax=float(r), nrings=5, multi=6)
        cp = PlanoConvex(z_0, curvature, 0, n_glass, n_air, thickness, aperture)
        focal_length_cp = cp.focal_point() - z_0
        f_cp = cp.focal_point()
        op_cp = OutputPlane(f_cp + thickness/(1e6))
        elements_cp = [cp.side1(), cp.side2(), op_cp]
        rb_cp.propagate_bundle(elements_cp)
        dx_cp.append((lamda * focal_length_cp) / (2 * r))
        rms_cp.append(rb_cp.rms())

    plt.plot(r_values, rms_pc, 'r', label='RMS Plano-Convex')
    plt.plot(r_values, dx_pc, 'g', label='Diffraction Scale Plano-Convex')

    plt.plot(r_values, rms_cp, 'b', label='RMS Convex-Plano')
    plt.plot(r_values, dx_cp, 'g', label='Diffraction Scale Convex-Plano')  

    plt.xlabel('Input Beam Radius (D)')
    plt.title('RMS and Diffraction Scale vs Input Beam Radius')
    plt.legend()

    rb_pc_35 = RayBundle(rmax=3.5, nrings=5, multi=6)
    pc_35 = PlanoConvex(z_0, 0, -1 * curvature, n_glass, n_air, thickness, aperture)
    focal_length_pc_35 = pc_35.focal_point() - z_0
    f_pc_35 = pc_35.focal_point()
    op_pc_35 = OutputPlane(f_pc_35)
    elements_pc_35 = [pc_35.side1(), pc_35.side2(), op_pc_35]
    rb_pc_35.propagate_bundle(elements_pc_35)
    rms_pc_35 = rb_pc_35.rms()
    dx_pc_35 = (lamda * (focal_length_pc_35 - thickness)) / (2 * 3.5)

    rb_cp_35 = RayBundle(rmax=3.5, nrings=5, multi=6)
    cp_35 = PlanoConvex(z_0, curvature, 0, n_glass, n_air, thickness, aperture)
    focal_length_cp_35 = cp_35.focal_point() - z_0
    f_cp_35 = cp_35.focal_point()
    op_cp_35 = OutputPlane(f_cp_35 + thickness/(1e6))
    elements_cp_35 = [cp_35.side1(), cp_35.side2(), op_cp_35]
    rb_cp_35.propagate_bundle(elements_cp_35)
    rms_cp_35 = rb_cp_35.rms()

    return fig, rms_pc_35, rms_cp_35, dx_pc_35

task16()

if __name__ == "__main__":

    # Run task 8 function
    task8()

    # Run task 10 function
    # FIG10 = task10()

    # Run task 11 function
    # FIG11, FOCAL_POINT = task11()

    # Run task 12 function
    # FIG12 = task12()

    # Run task 13 function
    # FIG13, TASK13_RMS = task13()

    # Run task 14 function
    # FIG14, TASK14_RMS, TASK14_DIFF_SCALE = task14()

    # Run task 15 function
    # FIG15_PC, FOCAL_POINT_PC, FIG15_CP, FOCAL_POINT_CP = task15()

    # Run task 16 function
    # FIG16, PC_RMS, CP_RMS, TASK16_DIFF_SCALE = task16()

    plt.show()
