import pickle
import sys

from abaqus import *
from abaqusConstants import C3D8H, DEFORMABLE_BODY, FIXED, INSTANTANEOUS, OFF, ON, PRONY, STEP, THREADS, THREE_D, TIME
from abaqusConstants import COMBINED, PARAMETERS

import regionToolset
import visualization
import mesh
import step                                               # noinspection UnusedImport

import numpy as np

package_directory = sys.argv[-1]
if package_directory not in sys.path:
    sys.path.append(package_directory)

from materials.SS2506 import ElasticPlasticTransformMaterial       # noqa

with open('parameter_pickle.pkl', 'r') as parameter_pickle:
    exx, eyy, ezz, pxx, pyy, pzz = pickle.load(parameter_pickle)
    material = pickle.load(parameter_pickle)            # type: ElasticPlasticTransformMaterial
    time_period = pickle.load(parameter_pickle)

backwardCompatibility.setValues(includeDeprecated=True, reportDeprecated=False)

# Cube 1 x 1 x 1
lx = 1.
ly = 1.
lz = 1.

# =============================================================================
#                          Creating the cube
# =============================================================================

model = mdb.Model(name='Cube')
sketch = model.ConstrainedSketch(name='Sketch S1', sheetSize=1)

p1 = (0.0, 0.0)
p2 = (lx, 0.0)
p3 = (lx, ly)
p4 = (0.0, ly)

sketch.Line(point1=p1, point2=p2)
sketch.Line(point1=p2, point2=p3)
sketch.Line(point1=p3, point2=p4)
sketch.Line(point1=p4, point2=p1)

part = model.Part(name='cube',
                  dimensionality=THREE_D,
                  type=DEFORMABLE_BODY)

part.BaseSolidExtrude(sketch=sketch, depth=lz)

# ===========================================================================
#                               Assembly
# ===========================================================================

assembly = model.rootAssembly
instance = assembly.Instance(name='Instance', part=part, dependent=ON)

part.seedEdgeBySize(edges=part.edges,
                    size=max(lx, ly, lz),
                    constraint=FIXED)
e_type = mesh.ElemType(elemCode=C3D8H)
part.setElementType(regions=(part.cells,),
                    elemTypes=(e_type,))
part.generateMesh()

# ==========================================================================
#                                Materials
# ==========================================================================

mat1 = model.Material(name='material')
mat1.Elastic(table=((material.E, material.v),),
             moduli=INSTANTANEOUS)

hardening_table = [material.sy0]
backstresses = len(material.gamma_m)
for i in range(backstresses):
    hardening_table += [material.Cm[i], material.gamma_m[i]]
mat1.Plastic(hardening=COMBINED, dataType=PARAMETERS, numBackstresses=backstresses,
             table=(hardening_table, ))
mat1.plastic.CyclicHardening(parameters=ON, table=((material.sy0, material.Q, material.b), ))

# ===========================================================================
#                                Sections
# ===========================================================================

model.HomogeneousSolidSection(name='Section',
                              material='material')
reg = (part.elements,)
part.SectionAssignment(region=reg,
                       sectionName='Section')

# ==========================================================================
#                                  Step
# ==========================================================================

comp = model.StaticStep(name='step',
                        previous='Initial',
                        timePeriod=time_period,
                        initialInc=1E-5*time_period,
                        maxNumInc=10000000,
                        minInc=1E-12*time_period,
                        maxInc=1E-2*time_period,
                        nlgeom=OFF)

comp.control.setValues(allowPropagation=OFF,
                       lineSearch=(5.0, 1.0, 0.0001, 0.25, 0.1),
                       resetDefaultValues=OFF)

# ==========================================================================
#                           Boundary conditions
# ==========================================================================

# Faces that cannot move, negative normal
xf = instance.faces.findAt(((0, 0.5*ly, 0.5*lz),),)
reg = regionToolset.Region(faces=xf)
model.DisplacementBC(name='dx1', createStepName='step', region=reg, u1=0)

yf = instance.faces.findAt(((0.5*lx, 0, 0.5*lz),),)
reg = regionToolset.Region(faces=yf)
model.DisplacementBC(name='dy1', createStepName='step', region=reg, u2=0)

zf = instance.faces.findAt(((0.5*lx, 0.5*ly, 0.),),)
reg = regionToolset.Region(faces=zf)
model.DisplacementBC(name='dz1', createStepName='step', region=reg, u3=0)

# Faces that can move, positive normal


def create_bc(direction, e, p):
    if direction is 'x':
        face_point = ((1.*lx, 0.5*ly, 0.5*lz),)
        bc_par = {'u1': lx if e is not None else None}
    elif direction is 'y':
        face_point = ((0.5*lx, 1.*ly, 0.5*lz),)
        bc_par = {'u2': ly if e is not None else None}
    elif direction is 'z':
        face_point = ((0.5*lx, 0.5*ly, 1.*lz),)
        bc_par = {'u3': lz if e is not None else None}
    else:
        raise (ValueError('Invalid axis'))

    face = instance.faces.findAt(face_point)
    region = regionToolset.Region(faces=face)
    surf = model.rootAssembly.Surface(name=direction + 'surf',
                                      side1Faces=face)

    if e is not None:
        model.TabularAmplitude(name='e' + direction + direction + '_amp',
                               timeSpan=STEP,
                               data=e.tolist())
        bc_par['name'] = 'd' + direction + '2'
        bc_par['createStepName'] = 'step'
        bc_par['region'] = region
        bc_par['amplitude'] = 'e' + direction + direction + '_amp'
        model.DisplacementBC(**bc_par)
    elif p is not None:
        model.TabularAmplitude(name='p' + direction + direction + '_amp',
                               timeSpan=STEP,
                               data=p.tolist())

        model.Pressure(name='p' + direction + '2',
                       createStepName='step',
                       region=surf,
                       magnitude=p,
                       amplitude='p' + direction + direction + '_amp')


create_bc('x', exx, pxx)
create_bc('y', eyy, pyy)
create_bc('z', ezz, pzz)

# ==========================================================================
#                           Output Requests
# ==========================================================================
model.FieldOutputRequest(name='F-Output-1',
                         createStepName='step',
                         frequency=1,
                         variables=('COORD',
                                    'U',
                                    'S',
                                    'E',
                                    'RF',))

# ==========================================================================
#                                 Job
# ==========================================================================

job = mdb.Job(name='oneElement',
              model=model,
              numCpus=1,
              numDomains=1,
              multiprocessingMode=THREADS)

job.submit()
job.waitForCompletion()

# ==========================================================================
#                           Output to file
# ==========================================================================

odb = visualization.openOdb('oneElement.odb')
frames_load = odb.steps['step'].frames

data = np.zeros((len(frames_load), 13))
for i in range(0, len(frames_load)):
    frame = frames_load[i]
    S = frame.fieldOutputs['S'].values[0]
    E = frame.fieldOutputs['E'].values[0]
    data[i, 0:6] = S.data
    data[i, 6:12] = E.data
    data[i, 12] = frame.frameValue

with open('stressStrain.pkl', 'wb') as pickle_handle:
    pickle.dump(data, pickle_handle)
odb.close()
