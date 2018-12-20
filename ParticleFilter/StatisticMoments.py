"""
Created on Fri Nov 03 09:40:51 2017

Author: Masha Itkina
Collaborator: Henry Shi

Compute the cell velocity statistics (means and variances).
"""

import numpy as np
import pdb

# INPUT:
# particle_array: object of class ParticleArray
def StatisticMoments(particle_array, grid_cell_array):

    state_size = particle_array.get_state_size()
    if state_size == 4:
        
		# arrays of length C
        vel_x_array_accum = particle_array.accumulate_cell_state(grid_cell_array, index_in_state = 2)
        vel_y_array_accum = particle_array.accumulate_cell_state(grid_cell_array, index_in_state = 3)
        vel_x_squared_array_accum = \
            particle_array.accumulate_cell_state_second_order(grid_cell_array, index_in_state1 = 2, index_in_state2 = 2)
        vel_y_squared_array_accum = \
            particle_array.accumulate_cell_state_second_order(grid_cell_array, index_in_state1 = 3, index_in_state2 = 3)
        vel_xy_array_accum = \
            particle_array.accumulate_cell_state_second_order(grid_cell_array, index_in_state1 = 2, index_in_state2 = 3)
    
        for j in range(grid_cell_array.get_length()):
            rho_p = grid_cell_array.get_cell_attr(j, "rho_p")
            start_idx = grid_cell_array.get_cell_attr(j, "start_index")
            end_idx = grid_cell_array.get_cell_attr(j, "end_index")

            if start_idx != None:
                # calc mean
                mean_x_vel = 1. / rho_p * vel_x_array_accum[j]
                mean_y_vel = 1. / rho_p * vel_y_array_accum[j]
                # calc variance
                
                var_x_vel = 1. / rho_p * (vel_x_squared_array_accum[j]) - mean_x_vel**2
                var_y_vel = 1. / rho_p * (vel_y_squared_array_accum[j]) - mean_y_vel**2
		
                # calc covariance
                covar_xy_vel = 1. / rho_p * (vel_xy_array_accum[j]) - mean_x_vel * mean_y_vel
                grid_cell_array.store_cell_stat_mean_vel_values(j, mean_x_vel, mean_y_vel, var_x_vel, var_y_vel, covar_xy_vel)

            else:
                grid_cell_array.store_cell_stat_mean_vel_values(j, 0., 0., 0., 0., 0.)
            
    else:
        raise Exception("Unexpected state size.")

    return
