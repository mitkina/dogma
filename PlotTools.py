#!/usr/bin/env python
"""Creates heading plots in the same style as Nuss et. al.

Nuss et. al. use a color wheel to represent the direction of
an object's velocity. This file creates velocity plots, then
overlays a circular colorbar legend.

The plotting tools also support displaying the environment,
such as unknown or occupied grid cells.

Author: Michael Anderson
"""

import matplotlib.pyplot as plt
import numpy as np
from matplotlib import cm, colors
import matplotlib as mpl

# Just for test cases
import sys
sys.path.insert(0, 'ParticleFilter/')
from Simulator import Grid, Car, Pedestrian


## Core Functionality
## DOGMAs
def colorwheel_plot(head_grid, occ_grid=None, m_occ_grid=None,
                    title="", show=True, save=False):
    """Create velocity plot of grid with color wheel overlay.
    
    INPUTS:
        head_grid - (np.matrix) Grid of the heading angles (rad)
        occ_grid - (opt)(np.matrix) Grid of occupancy statuses (0,1,2)
        title - (opt)(str) Title of generated plot
        show - (opt)(bool) Bool to display plot to user
        save - (opt)(bool) Bool to save plot to file
    OUTPUTS:
        None
    
    DESCRIPTION:
    head_grid is a matrix storing the *HEADING ANGLE* (rad) at each
    index. A heading angle of 0rad corresponds to moving right
    on the plot. This behavior can be controlled by the value
    of set_theta_offset().

    If an object has negligible velocity, you can opt not to plot
    any colors (instead of defaulting to 0rad heading). To force
    these indices to render as white, set their value to None.

    If you also want to see the environment plotted, we can overlay
    a second color map to display unknown/occluded/free space, stored
    in occ_grid. We assume the grid values are consistent with: 
        {unknown: 0, occupied: 1, free: 2}.
    """
    # Constants
    UNKNOWN_VAL = 0
    OCCUPIED_VAL = 1
    FREE_VAL = 2
    CAR = 3

    # Formatting
    lower_lim, upper_lim = [-np.pi, np.pi]
    color_map = 'hsv'
    # Frame percentages, formatted [left,bot,width,height]
    # TODO: Annoyingly, this is the percentage of full window (not the axes)
    legend_loc = [0.6,0.15,0.2,0.2]

    # Build the mask of unknown/occluded space
    # NOTE: unknown: 0, occupied: 1, free: 2
    if occ_grid is not None:
        env_cmap, env_norm = get_env_cmap()
        env_plt = plt.imshow(occ_grid, interpolation='none', \
                             cmap=env_cmap, norm=env_norm)

    # Build the mask of occupied masses
    # NOTE: occupied: 1, not occupied = 0
    #if m_occ_grid is not None:
        #occ_cmap = get_occ_cmap()
        #plt.clim(0, 2)
        #env_plt = plt.imshow(m_occ_grid, cmap=occ_cmap)

    # Create the original velocity plot
    # grid_plt = plt.imshow(grid, origin='lower')
    grid_plt = plt.imshow(head_grid, cmap=plt.cm.hsv, interpolation='none')
    plt.tick_params(axis='both', labelbottom='off', labelleft='off')
    plt.clim(lower_lim, upper_lim)
    if title:
        plt.title(title)

    # Overlay wheel legend
    legend_axes = plt.axes(legend_loc, projection='polar')
    legend_axes.set_theta_offset(np.pi/2.0) # Makes 0rad heading to the right
    # Wrap axes into circle; this hack may not be supported by some versions of Python
    legend_axes._direction = 2*np.pi  ## This is a nasty hack - using the hidden field to
                                      ## multiply the values such that 1 become 2*pi
                                      ## this field is supposed to take values 1 or -1 only!!
    norm = mpl.colors.Normalize(lower_lim, upper_lim)
    # Plot the colorbar onto the polar axis
    # note - use orientation horizontal so that the gradient goes around
    # the wheel rather than centre out
    quant_steps = 2056
    cb = mpl.colorbar.ColorbarBase(legend_axes,
                                   cmap=cm.get_cmap(color_map, quant_steps),
                                   norm=norm, orientation='horizontal')
    # aesthetics - get rid of border and axis labels  
    cb.outline.set_visible(False)
    legend_axes.set_axis_off()
    legend_axes.set_rlim([-1,1])
    plt.title("") # Prevent color wheel from also sharing title
    if save:
        if not title:  title = "plot"
        # else:  title = title.replace(" ", "_")
        fname = title + '.png'
        plt.savefig(fname)
    if show:
        plt.show()
    plt.close()
    return


def get_env_cmap():
    """Generate custom color map for UNKNOWN/OCCUPIED/FREE."""
    cmap = colors.ListedColormap(["#7575a3", "#000000", "#ffffff", "#000000"]) # grey, black, white, black
    bounds = [-0.5,0.5,1.5,2.5,3.5] # 0 = unknown (grey), 1 = occupied (black), 2 = free (white), 3 = car (black)
    norm = colors.BoundaryNorm(bounds, cmap.N)
    return cmap, norm

def get_occ_cmap():
    """Generate custom color map for UNKNOWN/OCCUPIED/FREE."""
    cmap = 'binary' # grey, black, white
    #norm = colors.BoundaryNorm(bounds, cmap.N)
    return cmap#, norm


## Particles
def particle_plot(particle_array, epsilon=1e-4,
                  title="", show=True, save=False, xlim = None, ylim = None):
    """Scatter/Quiver plot displaying particle positions/velocities.
    
    INPUTS:
        particle_array - (ParticleArray) container of Particles
        epsilon - (opt)(float) Minimum cell vel mag required to plot vel arrow
        title - (opt)(str) Title of generated plot
        show - (opt)(bool) Bool to display plot to user
        save - (opt)(bool) Bool to save plot to file
    OUTPUTS:
        None
    """
    assert epsilon > 0, "Quiver plots do not support velocities of zero."
    # Aesthetic parameters
    marker_dot_size = 9
    arrow_head_width = 3
    arrow_head_length = 5

    # Scatter plot of particle positions
    x_pos = np.array([particle[0] for particle in particle_array])
    y_pos = np.array([particle[1] for particle in particle_array])
    plt.plot(x_pos, y_pos, 'ko', markersize=marker_dot_size)
    plt.grid(True)

    # Quiver plot of particle velocities
    vel_x = np.array([particle[2] for particle in particle_array])
    vel_y = np.array([particle[3] for particle in particle_array])
    vel_mag = np.array([np.sqrt(vx**2+vy**2) for vx,vy in zip(vel_x,vel_y)])
    
    # Grab the indices with non-zero velocities magnitudes
    quiv_inds = np.where(np.array(vel_mag) > epsilon)
    if quiv_inds[0].size > 0:
        quiv_x, quiv_y = x_pos[quiv_inds], y_pos[quiv_inds]
        quiv_vel_x, quiv_vel_y = vel_x[quiv_inds], vel_y[quiv_inds]
        quiv_vel_mag = vel_mag[quiv_inds]
        quiv_heads = [np.arctan2(vy,vx) for vx,vy in zip(quiv_vel_x,quiv_vel_y)]
        
        U = [mag*np.cos(head) for mag,head in zip(quiv_vel_mag,quiv_heads)]
        V = [mag*np.sin(head) for mag,head in zip(quiv_vel_mag,quiv_heads)]        
        Q = plt.quiver(quiv_x, quiv_y, U, V, scale_units='xy', 
                       headwidth=arrow_head_width, headlength=arrow_head_length)
        mean_vel = int(np.mean(quiv_vel_mag))
        # X, Y, Key Value, Label Text
        plt.quiverkey(Q, 0.82, 0.93, mean_vel, '%dm/s' % mean_vel, labelpos='E', coordinates='figure')

    # Add some extra margin to plot axes
    margin = 0.2
    axes = plt.gca()
    if xlim:
        xmin, xmax = xlim
    else:
        xmin, xmax = axes.get_xlim()
    if ylim:
        ymin, ymax = ylim
    else:
        ymin, ymax = axes.get_ylim()

    xdelta = xmax - xmin
    ydelta = ymax - ymin
    plt.xlim([xmin-margin*xdelta, xmax+margin*xdelta])
    plt.ylim([ymin-margin*ydelta, ymax+margin*ydelta])

    # Plot options
    if title:
        plt.title(title)
    if save:
        if not title:  title = "plot"
        else:  title = title.replace(" ", "_")
        fname = title + '.png'
        plt.savefig(fname)
    if show:
        plt.show()
    plt.close()
    return


## Test Cases
# These functions make certain assumptions to plot; e.g. the ego
# vehicle is the first vehicle in the grid's vehicle list.
# For real implementation, these restrictions will be removed.

def sim_head_grid(grid, steps=10):
    """Generates heading grid from the simulator."""
    grid.tick(steps=steps, do_plot=False)
    head_grid = np.full(grid.global_grid.shape, None, dtype=float)
    for veh_num, vehicle in enumerate(grid.vehicles):
        if veh_num == 0:  # don't plot ego vehicle (no velocity assumed)
            continue
        # Compute heading and fill vehicle cells with heading value
        heading = vehicle.get_heading() # rad
        cells = grid.get_vehicle_cells(vehicle)
        for i,j in cells:
            head_grid[i,j] = heading
    return head_grid

    
def setup_simulator():
    """Creates mock instance of the simulator for testing."""
    # Define environment size (0.33 m resolution)
    grid_side_len = 256
    grid_spacing = 1

    # Define vehicles
    vehicles = []

    # Define stationary ego vehicle
    mid_pos = 1.0 * (grid_side_len / 2) / grid_spacing
    ego_pos = [mid_pos, mid_pos]
    ego_vel = [0.0,0.0]
    ego = Car(ego_pos, ego_vel)
    vehicles.append(ego)

    # Define other stuff
    car1 = Car([10.0,108.0], [0.0,1.0])
    vehicles.append(car1)

    ped1 = Pedestrian([150.0,148.0], [-1.0,-0.5])
    vehicles.append(ped1)

    ped2 = Pedestrian([50.0,50.0],[1.0,-0.5])
    vehicles.append(ped2)
    
    counter = 0
    grid = Grid(vehicles, grid_side_len, grid_spacing, counter)
    return grid


## Main

def main():
    grid = setup_simulator()
    steps_per_iter = 5
    for i in xrange(3):
        head_grid = sim_head_grid(grid, steps=steps_per_iter)

        title = "Simulation Iteration %d" % ((i+1) * steps_per_iter)
        #colorwheel_plot(vel_grid, title, show=False, save=True)
        colorwheel_plot(head_grid, occ_grid=grid.global_grid, title=title,
                        show=False, save=True)
    return

if __name__ == '__main__':
    seed = 1919 # Random.org
    np.random.seed(seed)
    main()
