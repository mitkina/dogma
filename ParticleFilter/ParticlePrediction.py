"""
Created on Fri Nov 03 09:40:51 2017

Author: Masha Itkina

Propagates the particles according to transition model (linear).
"""

from Particle import *
from Grid import *
import pdb

# INPUT:
# particle_array: class ParticleList object
# used to have grid_cell_array as input
def ParticlePrediction(particle_array, grid_cell_array, res):
    # get grid dimensions from grid_cell_array to map particles to the right index
    X, Y = grid_cell_array.get_shape()

    for i in range(particle_array.get_length()):
        # Applies equation 14
        particle_array.predict(i)

        # Replacing particle that goes off grid with the same particle in a uniformly random location in the grid
        particle_state = particle_array.get_state(i)
        x, y = particle_state[0:2]

        if (x > X - 1 or x < 0) or (y > Y - 1 or y < 0):
            particle_array.reinit(i)

        # grid cell needs to be between 0 and X - 1 or 0 and Y - 1
        particle_x = max(min(int(particle_state[0]), X - 1), 0)
        particle_y = max(min(int(particle_state[1]), Y - 1), 0)
        particle_array.set_cell_idx(i, particle_x * Y + particle_y) # used to be X

        # Applies equation 39
        particle_array.particles[i].persist_update_weight()
            
    return
