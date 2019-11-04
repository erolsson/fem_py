from materials.SS2506 import test_material


def write_input_file(filename, material, L=1, element_type='C3D8', umat_file=None, ):
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
                  '\t1, \t 0., 0., 0.',
                  '\t2, \t' + str(L) + ', 0., 0.',
                  '\t3, \t' + str(L) + ', ' + str(L) + ', 0.',
                  '\t4, \t 0., ' + str(L) + ', 0.',
                  '\t5, \t 0., 0., ' + str(L) + '',
                  '\t6, \t ' + str(L) + ', 0., ' + str(L) + '',
                  '\t7, \t ' + str(L) + ', ' + str(L) + ', ' + str(L) + '',
                  '\t8, \t 0., ' + str(L) + ', ' + str(L) + '',
                  '*Element, type=' + element_type + ', elset=all_elements',
                  '\t1, 1, 2, 3, 4, 5, 6, 7, 8',
                  '**',
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
        file_lines.append('\t\t' + material.umat_depvar)
        umat_parameters = material.umat_parameters
        file_lines.append('\t*User Material, constants=' + str(len(umat_parameters)))
        parameter_str = ''
        for i, par in enumerate(umat_parameters):
            parameter_str += str(par)
            if (i % 8 == 0 and i != 0) or i == len(umat_parameters) - 1:
                file_lines.append('\t\t' + parameter_str)
            else:
                parameter_str += ', '
    else:
        file_lines += material.abaqus_material_string

    with open(filename, 'w') as input_file:
        for line in file_lines:
            input_file.write(line + '\n')
        input_file.write('**EOF')


if __name__ == '__main__':
    write_input_file('test_inp_abq.inp', material=test_material)

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