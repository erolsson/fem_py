import pickle
import sys

import odbAccess

simulation_name = sys.argv[-1]

odb = odbAccess.openOdb(simulation_name + '.odb')
frames_load = odb.steps['step'].frames

data = np.zeros((len(frames_load), 13))
for i in range(0, len(frames_load)):
    frame = frames_load[i]
    S = frame.fieldOutputs['S'].values[0]
    E = frame.fieldOutputs['E'].values[0]
    data[i, 0:6] = S.data
    data[i, 6:12] = E.data
    data[i, 12] = frame.frameValue

with open('stressStrain' + simulation_name + '.pkl', 'wb') as pickle_handle:
    pickle.dump(data, pickle_handle)
odb.close()
