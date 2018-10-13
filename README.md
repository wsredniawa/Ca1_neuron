# Ca1_neuron
Python and hoc code to record and visualize neuron model and currents

Source of the Ca1 neuron model - https://senselab.med.yale.edu/ModelDB/ShowModel.cshtml?model=187610&file=/arrayTomography/#tabs-1

.txt files are example recordings of:
1. Rectime.txt - model integration time steps 
2. recstim.txt - current injectected to soma during one timestep (i.e. 10 nA)
3. matRecData.txt - all currents from the model: cap, pasive, Na, potassium, intracellular potential
4. matRecxyz.txt - position of neuron compartments in NEURON simulation
5. recsyn.txt - currents from synaptic mechanism. First row AMPA and second NMDA from all compartments.

Other files:
1. total_currents_ca1_CC.py - script to visualize recorded currents and cell morphology
2. save_hel.hoc - script that works after initialization of the original NEURON simulation and saves all currents from all cell segments
3. twinApical.swc - morphology of the cell


