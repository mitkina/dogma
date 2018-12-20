"""
Created on Fri Nov 03 09:40:51 2017

Author: Masha Itkina
Collaborator: Henry Shi

Resample particles from the newborn and persistent sets.
"""

import numpy as np
import copy
import pdb

def Resample(particle_array, birth_particle_array, particle_array_next, check_values = False):

    # check if birth particle array is empty
    if birth_particle_array.particles[0] == None:
        V = particle_array.get_length()
        rand_array = np.random.randint(V, size=V)

        # sort the array
        rand_array = np.sort(rand_array)

        # these weights should be normalized at the grid cell level (for each grid cell the weights should sum to rho_b or
        # rho_p respectively)
        persistent_accum = particle_array.accumulate_weight()
        offset = persistent_accum[-1]
        joint_weight_array_accum = persistent_accum

    else:
        V = particle_array.get_length()
        Vb = birth_particle_array.get_length()
        rand_array = np.random.randint(V + Vb, size = V)

        # sort the array
        rand_array = np.sort(rand_array)

        # sort the born particles array according to grid cell index
        birth_particle_array.sort_particles()

        # these weights should be normalized at the grid cell level (for each grid cell the weights should sum to rho_b or
        # rho_p respectively)
        persistent_accum = particle_array.accumulate_weight()
        offset = persistent_accum[-1]
        birth_accum = birth_particle_array.accumulate_weight() + offset
        joint_weight_array_accum =  np.concatenate((persistent_accum, birth_accum))

    # calculates resampled particle indices
    idx_array_resampled = calc_resampled_indeces(joint_weight_array_accum, rand_array, len(joint_weight_array_accum))

    # the weight for every particle must be the same and equal to sum of rho_b and rho_p for all particles divided by
    # number of particles
    new_particle_weight = joint_weight_array_accum[-1]/1./V

    for i in range(V):
        sample_index = idx_array_resampled[i]
        if sample_index < V:
            particle = particle_array.get_particle(sample_index)
        else:
            particle = birth_particle_array.get_particle(sample_index - V)

        # update to normalized equivalent weight
        particle.weight = new_particle_weight

        particle_array_next.set_particle(i, copy.deepcopy(particle))
        check_particle = particle_array_next.get_particle(i)

	# debugging information
    if check_values:
        weight_array_accum = particle_array_next.accumulate_weight()
        print "Final weight sum: ", weight_array_accum[-1]

        particle_array_next.sort_particles()

        # check weight array should be less than 1 in each grid cell
        index = particle_array_next.particles[0].index
        i = 0
        sum_weight = 0.
        while i != particle_array_next.get_length():
            if particle_array_next.particles[i].index == index:
                sum_weight += particle_array_next.particles[i].weight
                i = i + 1
            else:
                if (sum_weight <=1. and sum_weight >= 0.):
                    sum_weight = particle_array_next.particles[i].weight
                    index = particle_array_next.particles[i].index
                    i = i + 1

                else:
                    print "index: ", index, "particle: ", i
                    print "Sum of weights in cell: ", sum_weight
                    assert(sum_weight <= 1.)
                    assert(sum_weight >= 0.)

        print "Sum of weight in last set of particles: ", sum_weight

def calc_resampled_indeces(joint_weight_array_accum, rand_array, joint_particle_num):
    accum_max = joint_weight_array_accum[-1]
    scale_factor = (joint_particle_num - 1.) / accum_max
    weight_accum_scaled = joint_weight_array_accum * scale_factor

    # This should hold given the above calculation, but might not 
    # due to rounding error. We correct in case that happens.
    if weight_accum_scaled[-1] != rand_array[-1]: weight_accum_scaled[-1] = rand_array[-1]  
    
    # return the index of the range that incorporates the random selection
    random_particle_index = 0
    resampled_indeces = []
    for i in range(len(weight_accum_scaled)):
        while weight_accum_scaled[i] >= rand_array[random_particle_index]:
            resampled_indeces.append(i)
            random_particle_index += 1
            if random_particle_index == len(rand_array): 
                return resampled_indeces

    raise Exception("Should not get to here.")
