#!/usr/bin/env python
"""Simulates simple environment test case.
Author: Michael Anderson
"""

import numpy as np
import matplotlib.pyplot as plt
import pickle
import os

OUTPUT_DIR = '../data/sensor_grids/'

if not os.path.exists(OUTPUT_DIR):
	os.makedirs(OUTPUT_DIR)

def create_data(w=256, res=1):
    """Make an array for the demonstration."""
    n_elems = len(xrange(0, w, res))
    # unknown: 0, occupied: 1, free: 2
    Z = 2.0*np.ones((n_elems, n_elems))
    return Z


class Grid(object):
    """Represents GLOBAL environment grid."""
    def __init__(self, vehicles, w=256, dx=1, counter=0):  # FLAG
        self.vehicles = vehicles
        self.w, self.dx = w, dx
        self.counter = counter

        self.global_grid = create_data(self.w, self.dx)
        self.place_vehicles()

    def get_vehicle_cells(self, vehicle):
        """Return list of cells occupied by given vehicle."""
        center_x, center_y = vehicle.pos
        w, h = vehicle.dimens
        # Bounding box approximation
        left_edge = max(self.grid_floor(center_x - w/2.0), 0)
        right_edge = min(self.grid_ceil(center_x + w/2.0), self.global_grid.shape[0]-self.dx)
        bottom_edge = max(self.grid_floor(center_y - h/2.0), 0)
        top_edge = min(self.grid_ceil(center_y + h/2.0), self.global_grid.shape[1]-self.dx)
        # Map to grid indices
        edges = [left_edge, right_edge, bottom_edge, top_edge]
        l_ind, r_ind, b_ind, t_ind = [self.pos_to_ind(val) for val in edges]
        # Fill all encompassed indices
        vehicle_inds = []
        for i in xrange(l_ind, r_ind+1):
            for j in xrange(b_ind, t_ind+1):
                vehicle_inds.append((i,j))
        return vehicle_inds

    def grid_floor(self, val):
        """Floor to nearest grid value."""
        spacing = 1.0 / self.dx
        return int(np.floor(val * spacing)) / spacing

    def grid_ceil(self, val):
        """Ceil to nearest grid value."""
        spacing = 1.0 / self.dx
        return int(np.ceil(val * spacing)) / spacing

    def pos_to_ind(self, val):
        """Map precise position value to grid indix based on grid spacing."""
        ind = 1.0 * val / self.dx
        return int(ind)

    def place_vehicles(self):
        """Places vehicles onto global grid."""
        self.global_grid = create_data(self.w, self.dx)
        for vehicle in self.vehicles:
            inds = self.get_vehicle_cells(vehicle)
            for i,j in inds:
                self.global_grid[i,j] = 1
        return
    
    def create_global_grid(self, w=256, res=1.):
        """Creates x and y position global grids."""
        x_coords = np.arange(0,w*res,res)
        y_coords = np.arange(0,w*res,res)
        gridx,gridy = np.meshgrid(x_coords,y_coords)
      
        self.global_x_grid = gridx.T
        self.global_y_grid = gridy.T
        
        return

    def tick(self, steps=1, do_plot=False):
        for _ in xrange(steps):
            # Move all vehicles
            # .move() increments the *precise* internal position
            # .place_vehicles() assigns the vehicle to *coarser* grid indices
            [vehicle.move() for vehicle in self.vehicles]
            self.place_vehicles()

        if do_plot:
            plt.imshow(self.global_grid, \
                        cmap=plt.cm.binary, interpolation='none')
            plt.savefig(os.path.join(OUTPUT_DIR, 'simulation_' + str(self.counter) + '.png'))
            self.counter += 1
            plt.show()
        return
        


class Mover(object):
    """General interface for defining object dynamics."""
    def __init__(self, pos, vel, dimens):
        self.pos = np.array(pos)
        self.vel = np.array(vel)
        self.dimens = dimens
    def __repr__(self):
        return "Object at %s" % (self.pos,)

    # TODO: Add noise option
    def move(self, dt=0.1):
        self.pos += self.vel * dt

    def get_heading(self):
        return np.arctan2(self.vel[1], self.vel[0])

class Car(Mover):
    CAR_DIMENS = np.array([5,10]) # m
    def __init__(self, pos, vel):
        super(Car, self).__init__(pos, vel, Car.CAR_DIMENS)

    def __repr__(self):
        return "Car at %s" % (self.pos,)

class Pedestrian(Mover):
    PED_DIMENS = np.array([4,4]) # m
    def __init__(self, pos, vel):
        super(Pedestrian, self).__init__(pos, vel, Pedestrian.PED_DIMENS)

    def __repr__(self):
        return "Pedestrian at %s" % (self.pos,)


def main():
    # Define environment size (0.33 m resolution)
    grid_side_len = 256
    grid_spacing = 1

    # Define vehicles
    vehicles = []

    # Define stationary ego vehicle
    mid_pos = 1.0 * (grid_side_len / 2) / grid_spacing
    ego_pos = [mid_pos, mid_pos]
    ego_vel = [0,0]
    ego = Car(ego_pos, ego_vel)
    vehicles.append(ego)

    # Define other agents
    car1 = Car([100.25,100.], [0,30.])
    vehicles.append(car1)

    car2 = Car([150.25,150.], [0,-30.])
    vehicles.append(car2)
    
    # Define jpg labels counter
    counter = 0
    grids = []
    grid = Grid(vehicles, grid_side_len, grid_spacing, counter)
    for _ in xrange(10):
        grid.tick(steps=1, do_plot=True)
        
        # save the data
        grids.append(grid.global_grid)
    
    # generate x and y position global grids
    grid.create_global_grid(256, 1.)
    
    # pickle for now for testing purposes
    with open(os.path.join(OUTPUT_DIR, 'simulation.pickle'), 'wb') as f:
        pickle.dump([grids, grid.global_x_grid, grid.global_y_grid], f)

if __name__ == "__main__":
    seed = 1987
    np.random.seed(seed)
    main()
