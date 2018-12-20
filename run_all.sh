#! /bin/bash

cd ParticleFilter

# Simulate an occupancy grid environment
python Simulator.py

# Run Particle Filter
python ParticleFilterFunctions.py
