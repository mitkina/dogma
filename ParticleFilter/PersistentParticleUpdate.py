"""
Created on Fri Nov 03 09:40:51 2017

Author: Masha Itkina
Collaborator: Henry Shi

Updating particles that are persistent.
"""

import numpy as np
import pdb

# INPUT:
# particle_array: class ParticleList object
# used to have grid_cell_array as input
def PersistentParticleUpdate(particle_array, grid_cell_array, meas_cell_array, check_values=False):

    weight_array = []

    for i in range(particle_array.get_length()):
        weight_array.append(update_unnorm(particle_array.get_weight(i), particle_array.get_cell_idx(i), meas_cell_array))

    weight_array_accum = accumulate_weight(weight_array)

    # Equation 73
    weight_array_original_accum = particle_array.accumulate_weight()

    for j in range(grid_cell_array.get_length()):
        start_index = grid_cell_array.get_cell_attr(j, "start_index")
        end_index = grid_cell_array.get_cell_attr(j, "end_index")

        if start_index != None:
            # cells with particles predicted into
            # predicted occupied mass
            if start_index == 0:
                m_occ_accum = weight_array_accum[end_index]
            else:
                m_occ_accum = weight_array_accum[end_index] - weight_array_accum[start_index - 1]

            # look up persistent part of posterior occupancy mass
            rho_p = grid_cell_array.get_cell_attr(j, "rho_p")

            # normalization component for associated measurements
            mu_A = calc_norm_assoc(m_occ_accum, rho_p)

            if start_index == 0:
                m_occ_accum_orig = weight_array_original_accum[end_index]
            else:
                m_occ_accum_orig = weight_array_original_accum[end_index] - weight_array_original_accum[start_index - 1]
            mu_UA = calc_norm_unassoc(rho_p, m_occ_accum_orig)
            grid_cell_array.set_cell_nc(j, mu_A, mu_UA)

        else:
            # cells without particles predicted into
            next

    particle_array.normalize_weight(weight_array, grid_cell_array, meas_cell_array)

	# debugging information
    if check_values:
        # check weight array should equal rho_p in each grid cell?
        index = 100
        rho_p = grid_cell_array.get_cell_attr(index, "rho_p")
        print "rho_p: ", rho_p
        start_index = grid_cell_array.get_cell_attr(index, "start_index")
        end_index = grid_cell_array.get_cell_attr(index, "end_index")
        weight_array_accum = particle_array.accumulate_weight()
        if (start_index != None and end_index != 1):
            print "Sum of particles in the grid cell: ", weight_array_accum[end_index] - weight_array_accum[start_index-1]

        print "This is the total weight of particles after norm: ", weight_array_accum[-1]
        print "This is the total weight of particles before norm: ", weight_array_original_accum[-1]

    return

# equation 72
def calc_norm_assoc(m_accum, rho_p):
    # handle m_accum = 0 cases
    mu_A = rho_p / m_accum if m_accum > 0. else 0.
    return mu_A

# equation 73
def calc_norm_unassoc(rho_p, m_occ):
    # handle m_occ = 0 cases
    mu_UA = rho_p / m_occ if m_occ > 0. else 0.
    return mu_UA

# equation 69
def update_unnorm(weight, index, meas_cell_array):

    # GET Gk+1 distribution: w = g*w
    meas_m_free, meas_m_occ = meas_cell_array.get_measurement(index)

	# pseudoG  = 1
    pseudoG = meas_cell_array.get_pseudoG()
    weight *= pseudoG

    return weight

# accumulated weights for particles
def accumulate_weight(weight_array):

    accumulated_array = np.zeros([len(weight_array)])
    accumulated_array[0] = weight_array[0]

    for i in range(1, len(accumulated_array)):
        accumulated_array[i] = accumulated_array[i - 1] + weight_array[i]

    return accumulated_array
