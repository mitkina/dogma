# -*- coding: utf-8 -*-
"""
Created on Sat Nov 04 16:21:18 2017

Author: Masha Itkina
Collaborator: Henry Shi

Particle class.

"""

import numpy as np
import Grid
import pdb

class Particle(object):
    """Represents a single particle."""
    # INPUT: 
    # state size: 4 (pos, vel)
    # T: time interval for measurements (10 Hz default)
    # p_S: persistence probability
    def __init__(self, position, index, weight, state_size=4, p_S=0.99, scale_vel = 5., scale_acc = 2., \
                 process_pos = 0.1, process_vel = 0.5, process_acc = 0.8):
        
        self.state_size = state_size
        self.p_S = p_S
        # NEED TO INITIALIZE: UNIFORM FOR NOW, WITH STATIC ASSUMPTION AT START
        self.weight = weight
        self.index = index
        self.associated = True

        # velocity, acceleration variance initialization
        self.scale_vel = scale_vel
        self.scale_acc = scale_acc

        # position, velocity, acceleration process noise
        self.process_pos = process_pos
        self.process_vel = process_vel
        self.process_acc = process_acc

        x, y = position
        if self.state_size == 4:
            self.state = np.zeros([state_size])
            self.state[0] = x
            self.state[1] = y
            self.state[2] = np.random.normal(loc=0.0, scale=self.scale_vel, size=None)
            self.state[3] = np.random.normal(loc=0.0, scale=self.scale_vel, size=None)
        else:
            raise Exception("Unexpected state size.")
                
    def __getitem__(self, index):
        return self.state[index]

    def predict(self, transition_matrix, noiseSD = 0.1):
        # process noise should be constant for each time step: VARIANCE WILL PROB NEED TO DIFFER FOR POSITION, VELOCITY, ACCELERATION
        zeta = np.zeros(self.state_size)
        zeta[0:2] = np.random.normal(0, self.process_pos, 2)
        if self.state_size == 4:
            zeta[2:] = np.random.normal(0, self.process_vel, 2) # length 4, while state is position, velocity
        # Equations 14
        self.state = np.dot(transition_matrix, self.state) + zeta
        return
    
    def persist_update_weight(self):
        self.weight = self.weight*self.p_S
        if (self.weight < 0.):
            print "ATTENTION: persistent update has negative weights in particle prediction"
        return

    def normalize(self, updated_weight, p_A, mu_A, mu_UA):
        self.weight = p_A * updated_weight * mu_A + (1. - p_A) * self.weight * mu_UA
        return

    
class ParticleArray(object):
    """Represents a list of particles."""
    # INPUT: 
    # nu: number of consistent particles
    # state_size: number of states (4 or 6)
    def __init__(self, nu = 2*10**6, grid_shape = (256, 256), state_size = 4, T = 0.1, p_S = 0.99, scale_vel = 5., scale_acc = 2., \
                 process_pos = 0.1, process_vel = 0.5, process_acc = 0.8, birth = False, empty_array = False):
        self.state_size = state_size
        self.T = T
        self.p_S = p_S
        # initialize a list of size nu of particles
        self.particles = []
        self.nu = nu
        self.grid_shape = grid_shape

        # velocity, acceleration variance initialization
        self.scale_vel = scale_vel
        self.scale_acc = scale_acc

        # position, velocity, acceleration process noise
        self.process_pos = process_pos
        self.process_vel = process_vel
        self.process_acc = process_acc

        if self.state_size == 4:
            # pos, vel transitions
            self.transition_matrix = np.array([[1., 0., self.T, 0.], \
                                               [0., 1., 0., self.T], \
                                               [0., 0., 1., 0.], \
                                               [0., 0., 0., 1.]])
        else:
            raise Exception("Unexpected state size.")

        if not empty_array: 
            X, Y = grid_shape
            C = X * Y
            if birth == True:
                weight = 0.
            else:
                weight = 1./nu
            for i in range(nu):
                index = np.random.choice(C)

                # initializing particles at the center of cells
                x = index % X + 0.5 # operators were switched before
                y = index / X + 0.5
                self.particles.append(Particle(position = (x, y), index = index, weight = weight, state_size = self.state_size, p_S = self.p_S, \
                                               scale_vel = self.scale_vel, scale_acc = self.scale_acc, process_pos = self.process_pos, \
                                               process_vel = self.process_vel, process_acc = self.process_acc))

        else:
            self.particles = [None for i in range(nu)]

    def __repr__(self):
        return self.particles

    def __getitem__(self, index):
        return self.particles[index]

    def __setitem__(self, index, particle):
        self.particles[index] = particle
        return

    def __iter__(self):
        return self.particles.__iter__()

    def __len__(self):
        return len(self.particles)

    # return length
    def get_length(self):
        return len(self.particles)
    
    # return state size
    def get_state_size(self):
        return self.state_size

    # sort by grid cell index        
    def sort_particles(self):
        # sort the list in place by cell index
        self.particles.sort(key=lambda x: x.index)
        return

    # return cell index
    def get_cell_idx(self, particle_index):
        return self.particles[particle_index].index

    # set cell index of particular particle
    def set_cell_idx(self, particle_index, cell_index):
        self.particles[particle_index].index = cell_index
        return
    
    # set cell index and also set associated flad to True
    def set_cell_idx_A(self, particle_index, cell_index):
        self.set_cell_idx(particle_index, cell_index)
        self.particles[particle_index].associated = True
        return

    # set cell index and also set associated flad to False
    def set_cell_idx_UA(self, particle_index, cell_index):
        self.set_cell_idx(particle_index, cell_index)
        self.particles[particle_index].associated = False
        return

    # returns weight of particular particle
    def get_weight(self, particle_index):
        return self.particles[particle_index].weight

    # set weight of particular particle
    def set_weight(self, particle_index, weight):
        self.particles[particle_index].weight = weight
        return

    # returns state of particular particle
    def get_state(self, particle_index):
        return self.particles[particle_index].state

    # sets state of particular particle
    def set_state(self, particle_index, state):
        self.particles[particle_index].state = state
        return

    # returns particle object
    def get_particle(self, particle_index):
        return self.particles[particle_index]

    # sets particle object
    def set_particle(self, particle_index, new_particle):
        self.particles[particle_index] = new_particle
        return


    # calls predict on a particular particle
    def predict(self, particle_index):
        self.particles[particle_index].predict(transition_matrix = self.transition_matrix)

    def reinit(self, particle_index):

        X, Y = self.grid_shape
        index = np.random.choice(X*Y)
        # initializing particles at the center of cells
        x = index % X + 0.5 # operators were switched before
        y = index / X + 0.5
        self.particles[particle_index] = \
                Particle(position=(x, y), index=index, weight=self.particles[particle_index].weight, state_size=self.state_size, p_S=self.p_S, \
                                               scale_vel = self.scale_vel, scale_acc = self.scale_acc, process_pos = self.process_pos, \
                                               process_vel = self.process_vel, process_acc = self.process_acc)

    # OUTPUT: accumulated weights for particles
    def accumulate_weight(self):
        accumulated_array = np.zeros([len(self.particles)])
        accumulated_array[0] = self.particles[0].weight
    
        for i in range(1, len(accumulated_array)):
            if (self.particles[i].weight < 0.):
                print "NEGATIVE WEIGHT ISSUE IN ARRAY ACCUMULATION"
            accumulated_array[i] = accumulated_array[i-1] + self.particles[i].weight
    
        return accumulated_array
    
    # INPUT: 
    # grid_cell_array: accumulated array will be the same shape as the particle_array
    # index_in_state: the index of the state variable that we are accumulating
    # OUTPUT: accumulated states in each grid cell
    def accumulate_cell_state(self, grid_cell_array, index_in_state):
        accumulated_array = np.zeros([grid_cell_array.get_length()])

        # loop through the grid cells
        for i in range(grid_cell_array.get_length()):
            
            # cell accumulation is 0. by default
            cell_accumulation = 0.
            start_particle_index = grid_cell_array.get_cell_attr(i, "start_index")
            end_particle_index = grid_cell_array.get_cell_attr(i, "end_index")
            if start_particle_index != None:
                # if cell has particles assignes, loop through 
                # the particles in the sorted particle array
                for j in range(start_particle_index, end_particle_index + 1):
                    if self.get_weight(j) < 0.:
                        print "THIS IS AN ISSUE", self.get_weight(j)
                    cell_accumulation += self.get_state(j)[index_in_state] * self.get_weight(j)
            
            # cell_accumulation gets written into the array, so 
            # this is not an overall cummulative array
            accumulated_array[i] = cell_accumulation
    
        return accumulated_array
    
    # INPUT: 
    # grid_cell_array: accumulated array will be the same shape as the particle_array
    # index_in_state1: the index of the state variable that we are accumulating
    # index_in_state2: the index of the state variable that we are accumulating
    # OUTPUT: accumulated second order states in each grid cell
    def accumulate_cell_state_second_order(self, grid_cell_array, index_in_state1, index_in_state2):
		accumulated_array = np.zeros([grid_cell_array.get_length()])

		# loop through the grid cells
		for i in range(grid_cell_array.get_length()):

			# cell accumulation is 0. by default
			cell_accumulation = 0.
			start_particle_index = grid_cell_array.get_cell_attr(i, "start_index")
			end_particle_index = grid_cell_array.get_cell_attr(i, "end_index")
			if start_particle_index != None:
				# if cell has particles assignes, loop through 
				# the particles in the sorted particle array
				for j in range(start_particle_index, end_particle_index + 1):
					cell_accumulation += self.get_state(j)[index_in_state1] * self.get_state(j)[index_in_state2] * self.get_weight(j)

			# cell_accumulation gets written into the array, so 
			# this is not an overall cummulative array
			accumulated_array[i] = cell_accumulation
    
		return accumulated_array

    # equation 74: initialize associated new-born particles
    def initialize_newborn_A(self, particle_index, measurement, position = (0., 0.)):
        """
        Not sure how equation 74 works, using a dummy distribution
        for now. Will need to change to actual distribution which
        should be a function of measurement.
        """

        newstate = np.pad(position, (0, self.state_size - len(position)), 'constant')
        #variance = np.diag(1., 1., self.scale_vel, self.scale_vel, self.scale_acc, self.scale_acc)
        #newstate = np.random.multivariate_normal(mean, variance)

        if self.state_size == 4:
            newstate[2] =  np.random.normal(loc=0.0, scale=self.scale_vel, size=None)
            newstate[3] = np.random.normal(loc=0.0, scale=self.scale_vel, size=None)
        else:
            raise Exception("Unexpected state size.")

        self.set_state(particle_index, newstate)
        return

    # equation 76: initialize unassociated new-born particles
    def initialize_newborn_UA(self, particle_index, position = (0., 0.)):
        """
        Not sure how equation 76 works, using a dummy distribution
        for now. Will need to change to actual distribution.
        """
        state_size = len(self.get_state(particle_index))
        mean = np.pad(position, (0, state_size - len(position)), 'constant')
        variance = np.diag(np.ones(state_size) * 0.1)
        newstate = np.random.multivariate_normal(mean, variance)
        self.set_state(particle_index, newstate)
        return

    # equation 71: normalizes particle weights
    def normalize_weight(self, updated_weight, grid_cell_array, meas_cell_array):
        for i in range(self.get_length()):
            p_A = meas_cell_array.get_cell_attr(self.get_cell_idx(i), "p_A")
            mu_A = grid_cell_array.get_cell_attr(self.get_cell_idx(i), "mu_A")
            mu_UA = grid_cell_array.get_cell_attr(self.get_cell_idx(i), "mu_UA")
            self.particles[i].normalize(updated_weight[i], p_A, mu_A, mu_UA)

	# if the sum of weights in a grid cell exceeds one, normalize to persistent probability
    def normalize_p_S(self, start_index, end_index):

        particle_sum = 0.

        for i in range(start_index, end_index + 1):
            particle_sum += self.particles[i].weight

        for i in range(start_index, end_index + 1):
            self.particles[i].weight = self.particles[i].weight / particle_sum * self.p_S		 

