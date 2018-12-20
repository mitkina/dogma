# -*- coding: utf-8 -*-
"""
Created on Sat Nov 04 18:42:16 2017

Author: Masha Itkina
Collaborator: Henry Shi
"""

import numpy as np

class GridCell(object):
    """Represents a single grid cell."""
    def __init__(self, x, y, m_free = 0., m_occ = 0.):
        
        self.start_index = None
        self.end_index = None     
        
        # spatial coodinates
        self.x = x
        self.y = y
        
        # new-born part of posterior occupancy mass
        self.rho_b = 0.
        # remaining persistent part of posterior occupancy mass
        self.rho_p = 0.
        
        # INITIAL MASS GOES ALL (1) TO UNKNOWN {F,O}
        # posterior free mass
        self.m_free = m_free
        # posterior occupied mass
        self.m_occ = m_occ
        
        # normalization component for associated measurements
        self.mu_A = 1.
        # normalization component for unassociated measurements
        self.mu_UA = 0.

        # statistical moments: velocity
        self.mean_x_vel = 0.
        self.mean_y_vel = 0.
        self.var_x_vel = 0.
        self.var_y_vel = 0.
        self.covar_xy_vel = 0.

        # statistical moments: acceleration
        self.mean_x_acc = 0.
        self.mean_y_acc = 0.
        self.var_x_acc = 0.
        self.var_y_acc = 0.
        self.covar_xy_acc = 0.
    
    def store_values(self, rho_b, rho_p, m_free, m_occ):

        # new-born part of posterior occupancy mass        
        self.rho_b = rho_b
        # remaining persistent part of posterior occupancy mass
        self.rho_p = rho_p
        # posterior free mass
        self.m_free = m_free
        # posterior occupied mass
        self.m_occ = m_occ
        
        return

    def store_stat_mean_vel_values(self, mean_x_vel, mean_y_vel, var_x_vel, var_y_vel, covar_xy_vel):

        self.mean_x_vel = mean_x_vel
        self.mean_y_vel = mean_y_vel
        self.var_x_vel = var_x_vel
        self.var_y_vel = var_y_vel
        self.covar_xy_vel = covar_xy_vel
        
        return

    def store_stat_mean_acc_values(self, mean_x_acc, mean_y_acc, var_x_acc, var_y_acc, covar_xy_acc):

        self.mean_x_acc = mean_x_acc
        self.mean_y_acc = mean_y_acc
        self.var_x_acc = var_x_acc
        self.var_y_acc = var_y_acc
        self.covar_xy_acc = covar_xy_acc

        return

    # stores normalization components
    def set_normalization_components(self, mu_A, mu_UA):
        self.mu_A = mu_A
        self.mu_UA = mu_UA
        return


class CellArray(object):
    """Generic class for arrays"""
    # INPUT: 
    # C: number of grid cells in a (global?) grid
    def __init__(self, shape = (256, 256)):
        # initialize the grid cell array
        self.shape = shape
        self.C = shape[0] * shape[1]
        self.cells = []
        
        X, Y = shape
        for i in range(self.C):
            x = i / X
            y = i % X
            self.cells.append(GridCell(x,y))

    # returns length
    def get_length(self):
        return self.C

    # returns attribute of particular cell
    def get_cell_attr(self, cell_index, attribute):
        return getattr(self.cells[cell_index], attribute)

    # sets attribute of a particular cell
    def set_cell_attr(self, cell_index, attribute, value):
        if (len(self.cells) <= cell_index):        
            print len(self.cells), cell_index
        setattr(self.cells[cell_index], attribute, value)
        return

    # stores value in particular cell
    def store_cell_values(self, cell_index, rho_b, rho_p, m_free, m_occ):
        self.cells[cell_index].store_values(rho_b, rho_p, m_free, m_occ)
        return

    # returns shape
    def get_shape(self):
        return self.shape


class GridCellArray(CellArray):
    """Represents the (global?) grid array."""
    # INPUT: 
    # C: number of grid cells in a (global?) grid
    # p_A: TODO CHOOSE p_A probability
    def __init__(self, shape = (256, 256), p_A = 1.0):
        super(GridCellArray, self).__init__(shape)
        # initialize Grid cell Array
        self.p_A = p_A
        self.born_masses = []
        for i in range(self.C):
            # initialize the born masses
            self.born_masses.append(0.)
            
    # INPUT: string name of GridCell attribute
    # returns array with accumulated GridCell attributes,
    # returned array will be the same shape as the grid
    def accumulate(self, attribute):
        n = len(self.cells)
        accumulated_array = np.zeros([n])
        # loop through the grid cells
        accumulated_array[0] = self.get_cell_attr(0, attribute)
        for i in range(1, n):
            accumulated_array[i] = accumulated_array[i-1] + self.get_cell_attr(i, attribute)
        return accumulated_array

    # reinitialize the start and end indices to None for every cell
    def init_indices(self):
        for i in range(self.C):
            self.cells[i].start_index = None
            self.cells[i].end_index = None
        return

    # sets normalization components for particular cell
    def set_cell_nc(self, cell_index, mu_A, mu_UA):
        self.cells[cell_index].set_normalization_components(mu_A, mu_UA)
        return

    # retuns p_A
    def get_p_A(self):
        return self.p_A

    # stores mean velocity values for particular cell
    def store_cell_stat_mean_vel_values(self, cell_index, mean_x_vel, mean_y_vel, var_x_vel, var_y_vel, covar_xy_vel):
        self.cells[cell_index].store_stat_mean_vel_values(mean_x_vel, mean_y_vel, var_x_vel, var_y_vel, covar_xy_vel)

    # stores mean acceleration values for particular cell
    def store_cell_stat_mean_acc_values(self, cell_index, mean_x_acc, mean_y_acc, var_x_acc, var_y_acc, covar_xy_acc):
        self.cells[cell_index].store_stat_mean_acc_values(mean_x_acc, mean_y_acc, var_x_acc, var_y_acc, covar_xy_acc)
            
class MeasCellArray(CellArray):
    """Represents the (global?) measurement grid array."""
    def __init__(self, meas_free, meas_occ, shape = (256, 256), pseudoG = 1.):
        super(MeasCellArray, self).__init__(shape)
        # equation 69
        self.pseudoG = pseudoG
        for i in range(self.C):
            # initialize the association probabilities, default to 0 
            self.set_cell_attr(i, "p_A", 1.0)
            # For measurement cells m_unknown =  1. -  m_free - m_occ
            self.store_cell_values(i, 0.0, 0.0, m_free = meas_free[i], m_occ = meas_occ[i])

    # return measurement for particular cell
    def get_measurement(self, cell_index):
        return (self.cells[cell_index].m_free, self.cells[cell_index].m_occ)

    # returns pseudoG
    def get_pseudoG(self):
        return self.pseudoG



                        
