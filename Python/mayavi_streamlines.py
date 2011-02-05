#!/usr/bin/python
from numpy import *
from enthought.mayavi.mlab import *
 
N = 100 # the bigger the number, the more calculations your script has to make, but you save MayaVi's Runge Kutta method a lot of trouble and the total plot time is lowered
a = 0. # the lowest coordinate
b = 1. # the highest coordinate
dt = b / N; # the step we need to make to get to each position on our grid
 
q = [1., -1., 1., -1.] # the values of our four charges
qpos = [[0.56, 0.56, 0.50], # and their positions
        [0.26, 0.76, 0.50],
        [0.66, 0.16, 0.50],
        [0.66, 0.86, 0.50]]
 
x,y,z = mgrid[a:b:dt, a:b:dt, 0.:1.:0.5] # here is the trick - we want a 2d plot, but MayaVi works best in 3D space. Let's just keep the z-axis to a minimum (two steps)
Ex, Ey, Ez = mgrid[a:b:dt, a:b:dt, 0.:1.:0.5] # create a few arrays to store the vector field in, using meshgrid
 
for i in range(N): # iterate over all rows
    for j in range(N): # and all columns
        Ex[i,j] = 0.0 # set the value of each point to 0, initially
        Ey[i,j] = 0.0
        for num in range(len(q)): # for each charge, calculate the electric field it provides
            rs = ((x[i,j] - qpos[num][0])**2 + (y[i,j] - qpos[num][1])**2) # the distance from the point to the charge, squared
            r = sqrt(rs)
            q1x = q[num] * (x[i,j] - qpos[num][0]) / (r * rs) # the x-component of the field
            q1y = q[num] * (y[i,j] - qpos[num][1]) / (r * rs) # and the y-component
            # this is $\frac{q}{r^2} \hat{\mathbf r}$ on component form
            Ex[i,j] = q1x + Ex[i,j] # now, add this to the electric field in this point, together with the contribution from the other charges
            Ey[i,j] = q1y + Ey[i,j]
 
fig = figure(fgcolor=(0,0,0), bgcolor=(1,1,1)) # set the background and foreground of our figure
#obj = quiver3d(x, y, z, Ex, Ey, Ez, line_width=1) # uncomment this if you want a quiver plot
streams = list() # create a list to hold all our streamlines (or flows if you speak MayaVi)
 
for s in range(len(q)): # for each charge, create a streamline seed
    stream = flow(x,y,z,Ex, Ey, Ez, seed_scale=0.5, seed_resolution=1, seedtype='sphere') # the seed resolution is set to a minimum initially to avoid extra calculations
    stream.stream_tracer.initial_integration_step = 0.01 # the integration step for the runge kutta method
    stream.stream_tracer.maximum_propagation = 20.0 # the maximum length each step should reach - lowered to avoid messy output
    stream.stream_tracer.integration_direction = 'both' # integrate in both directions
    stream.seed.widget.center = qpos[s] # set the stream widget to the same position as the charge
    stream.seed.widget.radius = dt * 2 # and its radius a bit bigger than the grid size
    stream.seed.widget.theta_resolution = 30 # make the resolution high enough to give a fair number of lines
    stream.seed.widget.phi_resolution = 1 # but we are looking at a plane for now, so let's not have any resolution in the z-direction
    stream.seed.widget.enabled = False # hide the widget itself
    streams.append(stream) # and eventually, add the stream to our list for convenience
 
xlab = xlabel("x") # set the labels
ylab = ylabel("y")
showm = show() # show everything
axesm = axes() # add some axes
axesm.axes.y_axis_visibility = False # remove the z-axis (named y for some MayaVi reason)
fig.scene.z_plus_view() # let's look at it from the top!
fig.scene.parallel_projection = True # and we don't need any projection when looking at it in 2D
print "Done."