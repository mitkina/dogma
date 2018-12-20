"""
Created on Fri Nov 03 09:40:51 2017

Author: Henry Shi
Collaborator: Masha Itkina

New particle generation.
"""

import Particle
import Grid
import numpy as np
import pdb

def NewParticleInitialization(Vb, grid_cell_array, meas_cell_array, birth_particle_array, check_values = False):
    # accumulated mass for new born particles by cell
    particle_orders_array_accum = grid_cell_array.accumulate("rho_b")

    # check if birth probability is zero
    if particle_orders_array_accum[-1] == 0.:
        # no particles will be initialized
        return

    # normalize number of particles to Vb
    normalize_particle_orders(particle_orders_array_accum, Vb)

    for j in range(grid_cell_array.get_length()):
        # getting the number of new-born particles assigned to each cell
        # some cells may not get any
        start_idx = calc_start_idx(particle_orders_array_accum, j)
        end_idx = calc_end_idx(particle_orders_array_accum, j)
        num_new_particles = end_idx - start_idx + 1 if start_idx <= end_idx else 0
        p_A = meas_cell_array.get_cell_attr(j, "p_A")
        nu_A = calc_num_assoc(p_A, num_new_particles)
        nu_UA = num_new_particles - nu_A

        # associated and unassociated weights
        w_A = calc_weight_assoc(nu_A, p_A, grid_cell_array, j)
        w_UA = calc_weight_unassoc(nu_UA, p_A, grid_cell_array, j)

        # position of the middle point of a cell, this will be used in 
        # initializing the position of the associated particle
        X, Y = grid_cell_array.get_shape()
        cell_position = (j / X + 0.5, j % X + 0.5)

        # initializing associated particles
        for i in range(start_idx, start_idx + nu_A):
            birth_particle_array.set_cell_idx_A(i, j)
            birth_particle_array.set_weight(i, w_A)
            initialize_new_particle_A(birth_particle_array, i, meas_cell_array, cell_position)
        for i in range(start_idx + nu_A, end_idx + 1):
            birth_particle_array.set_cell_idx_UA(i, j)
            birth_particle_array.set_weight(i, w_UA)
            initialize_new_particle_UA(birth_particle_array, i, cell_position)

	# debugging values
    if check_values:
        # check weight array should equal rho_b in each grid cell?
        index = birth_particle_array.particles[0].index
        i = 0
        sum_weight = 0.
        while i != birth_particle_array.get_length():
            if birth_particle_array.particles[i].index == index:
                sum_weight += birth_particle_array.particles[i].weight
                i = i+1
            else:
                if (sum_weight - grid_cell_array.get_cell_attr(index, "rho_b"))**2 < 10.**-10:
                    sum_weight = birth_particle_array.particles[i].weight
                    index = birth_particle_array.particles[i].index
                    i = i+1

                else:
                    print "Sum of weights vs rho_b: ", sum_weight, rho_b
                    assert ((sum_weight - grid_cell_array.get_cell_attr(index, "rho_b"))**2 < 10.**-10)

        rho_b = grid_cell_array.get_cell_attr(index, "rho_b")
        print "rho_b: ", rho_b
        print "Sum of particles in the grid cell: ", sum_weight

    return 

# normalizes particle_orders_array_accum so that the aggregating
# particles assigned to each grid cell is Vb particles in total
def normalize_particle_orders(particle_orders_array_accum, Vb):
    # maximum is the last element in accumulative array
    array_max = particle_orders_array_accum[-1]
    if array_max <= 0: raise Exception("Accumulative array is empty or negative.")
    particle_orders_array_accum *=  Vb / (1.0 * array_max)
    return

# Calculates first index in birth_particle_array of cell j
def calc_start_idx(particle_orders_array_accum, cell_index):
    if cell_index == 0:
        return 0
    else:
        return int(particle_orders_array_accum[cell_index - 1])

# Calculates last index in birth_particle_array of cell j
def calc_end_idx(particle_orders_array_accum, cell_index):
    # end_idx would be start_idx - 1 if mass = 0 for cell
    return int(particle_orders_array_accum[cell_index]) - 1

# equantion 79 and 80: calculates number of new associated particles
def calc_num_assoc(p_A, num_new_particles):
    return int(p_A * num_new_particles)

# equation 75: calculates weight of an associated new particle
def calc_weight_assoc(nu_A, p_A, grid_cell_array, cell_index):
    return p_A * grid_cell_array.get_cell_attr(cell_index, "rho_b") / (1.0 * nu_A) if nu_A > 0 else 0.

# equation 77: calculates weight of an unassociated new particle
def calc_weight_unassoc(nu_UA, p_A, grid_cell_array, cell_index):
    return (1 - p_A) * grid_cell_array.get_cell_attr(cell_index, "rho_b") / (1.0 * nu_UA) if nu_UA > 0 else 0.

# Initializes associated new-born particles according to equations 74
def initialize_new_particle_A(particle_array, particle_index, meas_cell_array, position):
    cell_index = particle_array.get_cell_idx(particle_index)
    particle_array.initialize_newborn_A(particle_index, meas_cell_array.get_measurement(cell_index), position)
    return

# Initializes unassociated new-born particles according to equations 76
def initialize_new_particle_UA(particle_array, particle_index, position):
    particle_array.initialize_newborn_UA(particle_index, position)
    return
