"""
Created on Fri Nov 03 09:40:51 2017

Author: Masha Itkina
Collaborator: Henry Shi

Update Dempster-Shafer masses for all cells.
"""

import Particle
import Grid
import numpy as np
import pdb

def OccupancyPredictionUpdate(meas_cell_array, grid_cell_array, particle_array, p_B, alpha, check_values = False):

    # accumulate the weights for each grid cell from the particle list (same size as grid)
    weight_array_accum = particle_array.accumulate_weight()

    for j in range(grid_cell_array.get_length()):
        if grid_cell_array.get_cell_attr(j, "start_index") != None:
            # predicted occupied mass
            if grid_cell_array.get_cell_attr(j, "start_index") == 0:
                m_occ_pred = weight_array_accum[grid_cell_array.get_cell_attr(j, "end_index")]
            else:
                m_occ_pred = weight_array_accum[grid_cell_array.get_cell_attr(j, "end_index")] - \
                    weight_array_accum[grid_cell_array.get_cell_attr(j, "start_index") - 1]

            # truncate to pS probability
            if (m_occ_pred > particle_array.p_S):
                m_occ_pred = particle_array.p_S
                particle_array.normalize_p_S(grid_cell_array.get_cell_attr(j, "start_index"), grid_cell_array.get_cell_attr(j, "end_index")) 

            # predicted free mass
            m_free_pred = predict_free_mass(grid_cell_array.get_cell_attr(j, "m_free"), m_occ_pred, alpha)

            # combine measurement and prediction to form posterior occupied and free masses
            m_occ_up, m_free_up = update_of(m_occ_pred, m_free_pred, \
                meas_cell_array.get_cell_attr(j, "m_occ"), meas_cell_array.get_cell_attr(j, "m_free"))

			# debugging information
            if check_values and ((m_occ_up == 1.) and (m_occ_pred > 0.5)):
                print "check updates", m_occ_up
                print "check predictions", m_occ_pred
                print "check actual predictions", meas_cell_array.get_cell_attr(j, "m_occ")
			
			# debugging information
            if check_values and (m_occ_up > 1 or m_occ_up < 0):
                print "start index: ", grid_cell_array.get_cell_attr(j, "start_index"), "end index: ", \
                grid_cell_array.get_cell_attr(j, "end_index"), "mass_occ: ", m_occ_up, "mass_free: ", \
                m_free_up, "weight[start]: ", weight_array_accum[grid_cell_array.get_cell_attr(j, "start_index") - 1], \
                "weight[end]: ", weight_array_accum[grid_cell_array.get_cell_attr(j, "end_index")]
                assert(m_occ_up <= 1. and m_occ_up >= 0.)
                assert (m_free_up <= 1. and m_free_up >= 0.)
                assert(m_occ_up + m_free_up <= 1.)

            # compute new-born part of posterior occupancy mass
            rho_b = separate_newborn_part(m_occ_pred, m_occ_up, p_B)

            # compute remaining persistent part of posterior occupancy mass (equation 68)
            rho_p = m_occ_up - rho_b

            if check_values: assert(all([rho_b >= 0., rho_p >= 0., m_free_up >= 0., m_occ_up >= 0., m_occ_up <=1., m_free_up <= 1.]))
            
            grid_cell_array.store_cell_values(j, rho_b, rho_p, m_free_up, m_occ_up)
        
        else:
            # cells with no particles.
            next
    
    return

# equation 62
def predict_free_mass(cell_m_free, m_occ_pred, alpha):
    # limited by the mass addition to 1
    # information aging (same as in static grid)
    m_free_pred = min(alpha * cell_m_free, 1. - m_occ_pred)
    return m_free_pred

# equation 63: perform dst update
def update_of(m_occ_pred, m_free_pred, meas_m_occ, meas_m_free):
    
    # predicted unknown mass
    m_unknown_pred = 1. - m_occ_pred - m_free_pred
    
    # measurement masses: meas_m_free, meas_m_occ
    meas_cell_unknown = 1. - meas_m_free - meas_m_occ
    
    # implement DST rule of combination
    K = m_free_pred * meas_m_occ + m_occ_pred * meas_m_free
    
    m_occ_up = (m_occ_pred * meas_cell_unknown + m_unknown_pred * meas_m_occ + m_occ_pred * meas_m_occ) / (1. - K)
    m_free_up = (m_free_pred * meas_cell_unknown + m_unknown_pred * meas_m_free + m_free_pred * meas_m_free) / (1. - K)    
    
    return m_occ_up, m_free_up

# equation 67: compute the newborn part of the mass
def separate_newborn_part(m_occ_pred, m_occ_up, p_B):
    
    rho_b = (m_occ_up*p_B*(1. - m_occ_pred))/(m_occ_pred + p_B*(1. - m_occ_pred))
    
    return rho_b
