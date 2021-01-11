from pymatgen.core.lattice import Lattice, get_points_in_spheres
from pymatgen.core.structure import Structure
from pymatgen.vis.structure_vtk import StructureVis
import vtk
from pymatgen.electronic_structure.core import Magmom
from pymatgen.core.sites import PeriodicSite


class VisSpinStruc(StructureVis):
    def __init__(self):
        super().__init__(show_polyhedron=False)

    def set_structure(self, structure , tags=None): # set structure to class NonCollStructure ?
        super().set_structure(structure, reset_camera=True, to_unit_cell=True)
        for site in structure:
            for sp, occu in sorted(site.species.items()):
                # spin = site.specie.spin
                spin = Magmom(site.properties.get('magmom', getattr(sp, 'spin', 0)))
                coords = site.coords
                # if spin != 0:
                #     print(coords, site.specie.spin)
                if float(spin) != 0:
                    self.add_spin(site, spin, structure.lattice, show_text_dir=True)

    def add_spin(self, site: PeriodicSite, spin: Magmom, lattice: Lattice, show_text_dir=False):

        def generate_random_co_vec(normalizedX):
            import random

            # from vtk if normalizedX is not SAXIS
            # https://vtk.org/Wiki/VTK/Examples/Python/GeometricObjects/Display/OrientedArrow
            # normalizedX = [0 for i in range(3)]
            # math.Subtract(endpoint, site_coor, normalizedX) # normalizedX is like saxis
            # length = math.Norm(normalizedX)
            # math.Normalize(normalizedX)

            math = vtk.vtkMath()
            normalizedY = [0 for i in range(3)]
            normalizedZ = [0 for i in range(3)]
            # The Z axis is an arbitrary vector cross X
            arbitrary = [0 for i in range(3)]
            arbitrary[0] = random.uniform(-10, 10)
            arbitrary[1] = random.uniform(-10, 10)
            arbitrary[2] = random.uniform(-10, 10)
            math.Cross(normalizedX, arbitrary, normalizedZ)
            math.Normalize(normalizedZ)

            # The Y axis is Z cross X
            math.Cross(normalizedZ, normalizedX, normalizedY)

            return normalizedY, normalizedZ



        def change_range(moment, max):
            if moment != 0:
                return moment.get_moment() * (max / abs(float(moment)))
            elif moment == 0:
                return moment

        # site_mag_mom_ncl = change_range(site['species'][0]['properties']['spin'], max=normalize_magmom)
        print('')
        print(site.frac_coords, spin.saxis)


        # initialize arrow specifications
        arrow = vtk.vtkArrowSource()
        arrow.SetShaftRadius(0.1)
        arrow.SetShaftResolution(50)
        arrow.SetTipLength(0.5)
        arrow.SetTipRadius(0.2)
        arrow.SetTipResolution(50)

        # matrix for orientation (only X matters, the rest is randomly generate to fulfill matrix properties)
        normalizedX = list(spin.saxis)
        normalizedY, normalizedZ = generate_random_co_vec(normalizedX)
        # # Create the direction cosine matrix
        matrix = vtk.vtkMatrix4x4()
        matrix.Identity()
        for i in range(3):
            matrix.SetElement(i, 0, normalizedX[i])
            matrix.SetElement(i, 1, normalizedY[i])
            matrix.SetElement(i, 2, normalizedZ[i])

        # Apply the transforms
        transform = vtk.vtkTransform()
        transform.Translate(site.coords)
        transform.Concatenate(matrix)

        # in order to adapt the size the spin (depending on magnitude or relative terms ?)
        # transform.Scale(length, length, length)

        # also -> change color if spin is +/- or in plane

        mapper = vtk.vtkPolyDataMapper()
        actor = vtk.vtkActor()
        mapper.SetInputConnection(arrow.GetOutputPort())
        actor.SetUserMatrix(transform.GetMatrix())
        actor.SetMapper(mapper)
        self.ren.AddActor(actor)

        show_text_dir=True
        if show_text_dir:
            from math import gcd
            try:
                n1 = int(round(spin.saxis[0], 3)*1000)
                n2 = int(round(spin.saxis[1], 3) * 1000)
                n3 = int(round(spin.saxis[2], 3) * 1000)
                gcd2 = gcd(n1, n2)
                gcd3 = gcd(gcd2, n3)
                mil_bra_indices = (int(n1/gcd3), int(n2/gcd3), int(n3/gcd3))

            except:
                # text = str(np.round(float(site['species'][0]['properties']['spin']), decimals=decimal_magn)) + '\n' + str(
                #     np.round(site['species'][0]['properties']['spin'].get_00t_magmom_with_xyz_saxis().saxis,
                #              decimals=decimal_dir))
                print("miller bravais indices cannot be found")
                raise NotImplementedError

            import warnings
            warnings.warn("changed manual scalling of add_text")

            text = '  <' + ''.join(map(str, mil_bra_indices)) + '>'
            self.add_text(site.coords, text, scale=0.2)


if __name__ == "__main__":
    from github_non_collinear_spin.NonColSpinStructure import NonCollinearSpinStructure
    struc = Structure.from_dict({'@module': 'pymatgen.core.structure', '@class': 'Structure', 'charge': None, 'lattice': {'matrix': [[0.0, 5.07217361747331, 5.07217361747331], [5.07217361747331, 0.0, 5.07217361747331], [5.07217361747331, 5.07217361747331, 0.0]], 'a': 7.173136720541758, 'b': 7.173136720541758, 'c': 7.173136720541758, 'alpha': 60.00000000000001, 'beta': 60.00000000000001, 'gamma': 60.00000000000001, 'volume': 260.9830654620022}, 'sites': [{'species': [{'element': 'Dy', 'occu': 1}], 'abc': [0.625, 0.125, 0.125], 'xyz': [1.2680434043683275, 3.804130213104983, 3.804130213104983], 'label': 'Dy', 'properties': {}}, {'species': [{'element': 'Dy', 'occu': 1}], 'abc': [0.125, 0.625, 0.125], 'xyz': [3.804130213104983, 1.2680434043683275, 3.804130213104983], 'label': 'Dy', 'properties': {}}, {'species': [{'element': 'Dy', 'occu': 1}], 'abc': [0.125, 0.125, 0.625], 'xyz': [3.804130213104983, 3.804130213104983, 1.2680434043683275], 'label': 'Dy', 'properties': {}}, {'species': [{'element': 'Dy', 'occu': 1}], 'abc': [0.125, 0.125, 0.125], 'xyz': [1.2680434043683275, 1.2680434043683275, 1.2680434043683275], 'label': 'Dy', 'properties': {}}, {'species': [{'element': 'Ti', 'occu': 1}], 'abc': [0.125, 0.625, 0.625], 'xyz': [6.340217021841638, 3.804130213104983, 3.804130213104983], 'label': 'Ti', 'properties': {}}, {'species': [{'element': 'Ti', 'occu': 1}], 'abc': [0.625, 0.125, 0.625], 'xyz': [3.804130213104983, 6.340217021841638, 3.804130213104983], 'label': 'Ti', 'properties': {}}, {'species': [{'element': 'Ti', 'occu': 1}], 'abc': [0.625, 0.625, 0.125], 'xyz': [3.804130213104983, 3.804130213104983, 6.340217021841638], 'label': 'Ti', 'properties': {}}, {'species': [{'element': 'Ti', 'occu': 1}], 'abc': [0.625, 0.625, 0.625], 'xyz': [6.340217021841638, 6.340217021841638, 6.340217021841638], 'label': 'Ti', 'properties': {}}, {'species': [{'element': 'O', 'occu': 1}], 'abc': [0.7055778297240636, 0.7055778297240636, 0.2944221702759364], 'xyz': [5.07217361747331, 5.07217361747331, 7.157626506000941], 'label': 'O', 'properties': {}}, {'species': [{'element': 'O', 'occu': 1}], 'abc': [0.2944221702759364, 0.2944221702759364, 0.7055778297240636], 'xyz': [5.07217361747331, 5.07217361747331, 2.9867207289456785], 'label': 'O', 'properties': {}}, {'species': [{'element': 'O', 'occu': 1}], 'abc': [0.7055778297240636, 0.2944221702759364, 0.7055778297240636], 'xyz': [5.07217361747331, 7.157626506000941, 5.07217361747331], 'label': 'O', 'properties': {}}, {'species': [{'element': 'O', 'occu': 1}], 'abc': [0.2944221702759364, 0.7055778297240636, 0.2944221702759364], 'xyz': [5.07217361747331, 2.9867207289456785, 5.07217361747331], 'label': 'O', 'properties': {}}, {'species': [{'element': 'O', 'occu': 1}], 'abc': [0.2944221702759364, 0.7055778297240636, 0.7055778297240636], 'xyz': [7.157626506000941, 5.07217361747331, 5.07217361747331], 'label': 'O', 'properties': {}}, {'species': [{'element': 'O', 'occu': 1}], 'abc': [0.7055778297240636, 0.2944221702759364, 0.2944221702759364], 'xyz': [2.9867207289456785, 5.07217361747331, 5.07217361747331], 'label': 'O', 'properties': {}}, {'species': [{'element': 'O', 'occu': 1}], 'abc': [0.5444221702759364, 0.5444221702759364, 0.9555778297240636], 'xyz': [7.608260426209965, 7.608260426209965, 5.522807537682334], 'label': 'O', 'properties': {}}, {'species': [{'element': 'O', 'occu': 1}], 'abc': [0.9555778297240636, 0.9555778297240636, 0.5444221702759364], 'xyz': [7.608260426209965, 7.608260426209965, 9.693713314737597], 'label': 'O', 'properties': {}}, {'species': [{'element': 'O', 'occu': 1}], 'abc': [0.5444221702759364, 0.9555778297240636, 0.5444221702759364], 'xyz': [7.608260426209965, 5.522807537682334, 7.608260426209965], 'label': 'O', 'properties': {}}, {'species': [{'element': 'O', 'occu': 1}], 'abc': [0.9555778297240636, 0.5444221702759364, 0.9555778297240636], 'xyz': [7.608260426209965, 9.693713314737597, 7.608260426209965], 'label': 'O', 'properties': {}}, {'species': [{'element': 'O', 'occu': 1}], 'abc': [0.9555778297240636, 0.5444221702759364, 0.5444221702759364], 'xyz': [5.522807537682334, 7.608260426209965, 7.608260426209965], 'label': 'O', 'properties': {}}, {'species': [{'element': 'O', 'occu': 1}], 'abc': [0.5444221702759364, 0.9555778297240636, 0.9555778297240636], 'xyz': [9.693713314737597, 7.608260426209965, 7.608260426209965], 'label': 'O', 'properties': {}}, {'species': [{'element': 'O', 'occu': 1}], 'abc': [0.25, 0.25, 0.25], 'xyz': [2.536086808736655, 2.536086808736655, 2.536086808736655], 'label': 'O', 'properties': {}}, {'species': [{'element': 'O', 'occu': 1}], 'abc': [0.0, 0.0, 0.0], 'xyz': [0.0, 0.0, 0.0], 'label': 'O', 'properties': {}}]})
    array = []
    for i in range(22):
        array.append((1, 1, 0))



    spined_structure = NonCollinearSpinStructure.from_sites_and_orientations(struc, array, spins_in_sites_or_species="sites")
    print(spined_structure)
    vis = VisSpinStruc()
    vis.set_structure(spined_structure)
    vis.show()

