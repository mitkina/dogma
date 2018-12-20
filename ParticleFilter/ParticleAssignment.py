"""
Created on Fri Nov 03 09:40:51 2017

Author: Masha Itkina
Collaborator: Henry Shi

Assigns particles to grid cell array.
"""

# inputs include weight_array
def ParticleAssignment(particle_array, grid_cell_array):
    
    # sort by grid_cell_index
    particle_array.sort_particles()

    # reinitialize the grid_cell_array start and end indices to be None
    grid_cell_array.init_indices()
    
    # set temp variables for keeping track of start
    prev_cell_idx = None
    
    for i in range(particle_array.get_length()):

        j = particle_array.get_cell_idx(i)
        
        # first particle
        if prev_cell_idx == None:
            grid_cell_array.set_cell_attr(j, "start_index", i)
            prev_cell_idx = j
        # the cell index changed
        elif prev_cell_idx != j:
            grid_cell_array.set_cell_attr(prev_cell_idx, "end_index", i-1)
            grid_cell_array.set_cell_attr(j, "start_index", i)
            prev_cell_idx = j
        else:
            next
        
        # sets end_index if at end of particle array
        if i == particle_array.get_length() - 1:
            grid_cell_array.set_cell_attr(j, "end_index", i)

		# debugging
        if particle_array.get_weight(i) <= 0:
            raise Exception("WHY ARE PARTICLE WEIGHTS ZERO")

    return
