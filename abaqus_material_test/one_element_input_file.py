from collections import namedtuple

import numpy as np

from materials.SS2506 import test_material

BC = namedtuple('BC', ['amplitude', 'direction', 'mode'])


def write_input_file(filename, material, boundary_conditions, element_size=1., element_type='C3D8', umat_file=None,
                     time_period=1., max_increment=1.):
    try:
        material_name = material.name
    except AttributeError:
        material_name = 'material'
    file_lines = ['**',
                  '*Heading',
                  '\t One element simulation for testing materials',
                  '*Preprint, echo=NO, model=NO, history=NO, contact=NO',
                  '**',
                  '** ----------------------------------------------------------------',
                  '**   Create Geometry',
                  '*Node, nset=all_nodes',
                  '\t1,\t0.,\t0.,\t0.',
                  '\t2,\t' + str(element_size) + ',\t0.,\t0.',
                  '\t3,\t' + str(element_size) + ',\t' + str(element_size) + ',\t0.',
                  '\t4,\t0.,\t' + str(element_size) + ',\t0.',
                  '\t5,\t0.,\t0.,\t' + str(element_size),
                  '\t6,\t' + str(element_size) + ',\t0.,\t' + str(element_size),
                  '\t7,\t' + str(element_size) + ',\t' + str(element_size) + ',\t' + str(element_size),
                  '\t8,\t0.,\t' + str(element_size) + ',\t' + str(element_size),
                  '*Element, type=' + element_type + ', elset=all_elements',
                  '\t1, 1, 2, 3, 4, 5, 6, 7, 8',
                  '**',
                  '*Surface, type=ELEMENT, name=xsurf',
                  '\t1, S6',
                  '*Surface, type=ELEMENT, name=ysurf',
                  '\t1, S5',
                  '*Surface, type=ELEMENT, name=zsurf',
                  '\t1, S2',
                  '**',
                  '*Nset, nset=xnodes',
                  '\t2, 3, 6, 7',
                  '*Nset, nset=ynodes',
                  '\t3, 4, 7, 8',
                  '*Nset, nset=znodes',
                  '\t5, 6, 7, 8',
                  '** ----------------------------------------------------------------',
                  '**',
                  '**   Define material properties',
                  '**',
                  '*Solid Section, elset=all_Elements, material=' + material_name,
                  '\t1.0',
                  '**',
                  '** DEFINE MATERIAL PROPERTIES',
                  '**',
                  '*Material, name=' + material_name]
    if umat_file:
        file_lines.append('\t*Depvar')
        file_lines.append('\t\t' + str(material.umat_depvar()))
        umat_parameters = material.umat_parameters()
        file_lines.append('\t*User Material, constants=' + str(len(umat_parameters)))
        parameter_str = ''
        for i, par in enumerate(umat_parameters):
            parameter_str += str(par)
            if (i % 8 == 0 and i != 0) or i == len(umat_parameters) - 1:
                file_lines.append('\t\t' + parameter_str)
                parameter_str = ''
            else:
                parameter_str += ', '
    else:
        file_lines += material.abaqus_material_string()

    # Fixed boundary conditions
    bc_lines = ['*Boundary',
                '\t1, 1, 3',
                '\t2, 2, 3',
                '\t3, 3, 3',
                '\t4, 1, 1',
                '\t4, 3, 3',
                '\t5, 1, 2',
                '\t6, 2, 2',
                '\t8, 1, 1',
                '**']
    file_lines += bc_lines

    # Amplitude definitions for the time dependent bc
    for bc in boundary_conditions:
        file_lines.append('*Amplitude, name=' + bc.direction + '_amp')
        for t, val in bc.amplitude:
            file_lines.append('\t' + str(t) + ', ' + str(val))

    file_lines.append('*Step, name=step, nlgeom=NO, inc=10000000')
    file_lines.append('\t*Static')
    file_lines.append('\t\t1e-05, ' + str(time_period) + ', 1e-12,  ' + str(max_increment))

    # Defining the boundary conditions
    dir_dict = {'x': '1', 'y': '2', 'z': '3'}
    for bc in boundary_conditions:
        if bc.mode == 'strain':
            file_lines.append('\t*Boundary, amplitude=' + bc.direction + '_amp')
            file_lines.append('\t\t' + bc.direction + 'nodes, ' + dir_dict[bc.direction] + ', ' + bc.direction + '1.')
        elif bc.mode == 'stress':
            file_lines.append('\t*Dsload, amplitude=' + bc.direction + '_amp')
            file_lines.append('\t\t' + bc.direction + 'surf, P, -1')
        else:
            raise ValueError("Mode for a boundary conditions must be either stress or strain")

    file_lines.append('\t*Output, field')
    file_lines.append('\t\t*Node Output')
    file_lines.append('\t\t\tCOORD, RF, U')
    file_lines.append('\t\t*Element Output, directions=YES')
    element_output = '\t\t\tE, S'
    if umat_file:
        element_output += ', SDV'
    file_lines.append(element_output)
    file_lines.append('*End Step')

    with open(filename, 'w') as input_file:
        for line in file_lines:
            input_file.write(line + '\n')
        input_file.write('**EOF')


if __name__ == '__main__':
    n = 100
    amp_data = np.zeros((n, 2))
    amp_data[:, 0] = np.linspace(0, 1., 100)
    amp_data[:, 1] = np.sin(2*np.pi*amp_data[:, 0])*3000
    bc_list = [BC(amplitude=amp_data, direction='x', mode='stress')]
    write_input_file('test_inp_abq.inp', material=test_material, boundary_conditions=bc_list)
    write_input_file('test_inp_umat.inp', material=test_material, umat_file=1, boundary_conditions=bc_list)

"""
*User Material, constants=15
200000.,    0.3,  1000.,  1000.,     0.,   100.,     3.,135000.
950.,700000.,   500., 50000.,    50.,     0.,     0.


** Section: Section
*Solid Section, elset=_M4, material=material
,
*End Part
**
**
** ASSEMBLY
**
*Assembly, name=Assembly
**
*Instance, name=Instance, part=cube
*End Instance
**
*Nset, nset=_PickedSet3, internal, instance=Instance, generate
2,  8,  2
*Elset, elset=_PickedSet3, internal, instance=Instance
1,
*Nset, nset=_PickedSet4, internal, instance=Instance, generate
1,  4,  1
*Elset, elset=_PickedSet4, internal, instance=Instance
1,
*Nset, nset=_PickedSet5, internal, instance=Instance
3, 4, 7, 8
*Elset, elset=_PickedSet5, internal, instance=Instance
1,
*Elset, elset=_xsurf_S6, internal, instance=Instance
1,
*Surface, type=ELEMENT, name=xsurf
_xsurf_S6, S6
*Elset, elset=_ysurf_S1, internal, instance=Instance
1,
*Surface, type=ELEMENT, name=ysurf
_ysurf_S1, S1
*Elset, elset=_zsurf_S3, internal, instance=Instance
1,
*Surface, type=ELEMENT, name=zsurf
_zsurf_S3, S3
*End Assembly
"""