import os
import pickle
import sys

import numpy as np

import odbAccess

simulation_directory = sys.argv[-2]
simulation_name = sys.argv[-1]

os.chdir(simulation_directory)
odb = odbAccess.openOdb(simulation_name + '.odb')
frames_load = odb.steps['step'].frames

data = np.zeros((len(frames_load), 15))
for i in range(0, len(frames_load)):
    frame = frames_load[i]
    S = frame.fieldOutputs['S'].values[0]
    E = frame.fieldOutputs['E'].values[0]
    epl = frame.fieldOutputs['SDV1'].values[0]
    fM = frame.fieldOutputs['SDV2'].values[0]
    data[i, 0:6] = S.data
    data[i, 6:12] = E.data
    data[i, 12] = epl.data
    data[i, 13] = fM.data
    data[i, 14] = frame.frameValue

with open('stressStrain' + simulation_name + '.pkl', 'wb') as pickle_handle:
    pickle.dump(data, pickle_handle)
odb.close()
