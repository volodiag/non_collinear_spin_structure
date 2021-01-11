from typing import Dict, List, Tuple, Optional, Union, Iterator, Set, Sequence, Iterable
from pymatgen.core.lattice import Lattice, get_points_in_spheres
from pymatgen.core.periodic_table import Element, Specie, get_el_sp, DummySpecie
from pymatgen.core.composition import Composition
from pymatgen.core.structure import Structure
from pymatgen.symmetry.analyzer import SpacegroupAnalyzer
# from pymatgen.vis.structure_vtk import StructureVis
# import vtk
from pymatgen.electronic_structure.core import Magmom
from pymatgen.core.periodic_table import Element
from pymatgen.core.periodic_table import Specie
import numpy as np
from pymatgen.core.sites import PeriodicSite
import os
from monty.serialization import loadfn, dumpfn
import pymatgen.analysis.magnetism.analyzer
import numpy as np
import warnings
from math import gcd
from pymatgen.io.cif import CifWriter
from collections import OrderedDict, deque
from github_non_collinear_spin.VisualizeNonColSpinStructure import VisSpinStruc

MODULE_DIR = os.path.dirname(os.path.abspath(pymatgen.analysis.magnetism.analyzer.__file__))
DEFAULT_MAGMOMS = loadfn(os.path.join(MODULE_DIR, "default_magmoms.yaml"))


def saxis_normalize_to_gcd(spin):
    print(spin)
    if float(spin) != 0:
        try:
            n1 = int(round(spin.saxis[0], 3) * 1000)
            n2 = int(round(spin.saxis[1], 3) * 1000)
            n3 = int(round(spin.saxis[2], 3) * 1000)
            gcd2 = gcd(n1, n2)
            gcd3 = gcd(gcd2, n3)
            mil_bra_indices = (int(n1 / gcd3), int(n2 / gcd3), int(n3 / gcd3))

        except:
            warnings.warn("miller bravais indices cannot be found")
            mil_bra_indices = spin.saxis
    else:
        mil_bra_indices = (0, 0, 0)
    return mil_bra_indices



class NonCollinearSpinStructure(Structure):
    def __init__(self,
                 lattice: Union[List, np.ndarray, Lattice],
                 species: Sequence[Union[str, Element, Specie, DummySpecie, Composition]],
                 ncl_spins: List[Magmom],
                 coords: Sequence[Sequence[float]],
                 charge: float = None,
                 validate_proximity: bool = False,
                 to_unit_cell: bool = True,
                 coords_are_cartesian: bool = False,
                 site_properties: dict = None,
                 auto_ox_state: bool = True,
                 spins_in_sites_or_species="sites",
                 ):
        super().__init__(
            lattice, species, coords, charge=charge,
            validate_proximity=validate_proximity, to_unit_cell=to_unit_cell,
            coords_are_cartesian=coords_are_cartesian,
            site_properties=site_properties)

        self._sites = list(self._sites)  # type: ignore


        if len(ncl_spins) != len(self._sites):
            raise ValueError("Spin of all sites must be "
                             "specified in the dictionary.")

        self.spins_type = spins_in_sites_or_species
        # non collinear spin is rarely shought with disorder -> species-based spin is tricky to implement and not further developped
        if spins_in_sites_or_species == "species":
            for site, spin in zip(self.sites, ncl_spins):
                new_sp = {}
                for sp, occu in site.species.items():
                    sym = sp.symbol
                    oxi_state = getattr(sp, "oxi_state", None)
                    new_sp[Specie(sym, oxidation_state=oxi_state,
                                  properties={'spin': spin})] = occu
                site.species = Composition(new_sp)
        else:
            if spins_in_sites_or_species == "sites":
                for site, spin in zip(self._sites, ncl_spins):
                    site.properties['magmom'] = spin
                    site.properties['orientation (for info)'] = '<' + ''.join(map(str, saxis_normalize_to_gcd(spin))) + '>'


    @classmethod
    def from_sites_and_orientations(cls, sites: List[PeriodicSite], #i.e. Structure instance
                                        orientations,  # list of (x, y, z) for each site
                                        charge: float = None, validate_proximity: bool = False,
                                        to_unit_cell: bool = False, user_magmom_by_el={}, spins_in_sites_or_species="site"):

        # same as "from_sites" until spin part
        if len(sites) < 1:
            raise ValueError("You need at least one site to construct a %s" %
                             cls)
        prop_keys = []  # type: List[str]
        props = {}
        lattice = sites[0].lattice
        for i, site in enumerate(sites):
            if site.lattice != lattice:
                raise ValueError("Sites must belong to the same lattice")
            for k, v in site.properties.items():
                if k not in prop_keys:
                    prop_keys.append(k)
                    props[k] = [None] * len(sites)
                props[k][i] = v
        for k, v in props.items():
            if any((vv is None for vv in v)):
                warnings.warn("Not all sites have property %s. Missing values "
                              "are set to None." % k)

        # spin part
        # enables magmom intensity twisting -> default magmom
        def get_magmom_int(user_magmom_by_el, site):
            if not user_magmom_by_el:
                try:
                    magmom = DEFAULT_MAGMOMS[str(list(site._species._data.keys())[0])]
                except KeyError:
                    magmom = 0
            else:
                raise NotImplemented
                # magmom = user_magmom.get(str(list(site._species._data.keys())[0]))
            return magmom

        #  construct list of magmom based on specie type and orientation
        #  orientations should belong to the same site -> is it possible to set via an interface ?
        ncl_spins = []
        for site, orientation in zip(sites, orientations):
            print(site, orientation)
            magmom_int = get_magmom_int(user_magmom_by_el, site) # function is right above
            magmom = Magmom((0, 0, magmom_int), orientation)
            # magmom = magmom.get_xyz_magmom_with_001_saxis()
            ncl_spins.append(magmom)

        return cls(lattice,
                   [site.species for site in sites],
                   ncl_spins,
                   [site.frac_coords for site in sites],
                   charge=charge,
                   site_properties=props,
                   validate_proximity=validate_proximity,
                   to_unit_cell=to_unit_cell,
                   spins_in_sites_or_species=spins_in_sites_or_species,)

    def set_spin_by_number(self, site_number, new_orientation : Tuple=False, new_intensity: float=False):
        if self.spins_type != "site":
            raise NotImplemented
        if new_orientation and new_intensity:
            self[site_number].properties['magmom'] = Magmom(moment=(0, 0, new_intensity), saxis=new_orientation).get_moment()
        else:
            warnings.warn("should precise both orientation and intensity")



def structure_pyrochlore_non_magn():
    return Structure.from_dict({'@module': 'pymatgen.core.structure', '@class': 'Structure', 'charge': None, 'lattice': {'matrix': [[0.0, 5.07217361747331, 5.07217361747331], [5.07217361747331, 0.0, 5.07217361747331], [5.07217361747331, 5.07217361747331, 0.0]], 'a': 7.173136720541758, 'b': 7.173136720541758, 'c': 7.173136720541758, 'alpha': 60.00000000000001, 'beta': 60.00000000000001, 'gamma': 60.00000000000001, 'volume': 260.9830654620022}, 'sites': [{'species': [{'element': 'Dy', 'occu': 1}], 'abc': [0.625, 0.125, 0.125], 'xyz': [1.2680434043683275, 3.804130213104983, 3.804130213104983], 'label': 'Dy', 'properties': {}}, {'species': [{'element': 'Dy', 'occu': 1}], 'abc': [0.125, 0.625, 0.125], 'xyz': [3.804130213104983, 1.2680434043683275, 3.804130213104983], 'label': 'Dy', 'properties': {}}, {'species': [{'element': 'Dy', 'occu': 1}], 'abc': [0.125, 0.125, 0.625], 'xyz': [3.804130213104983, 3.804130213104983, 1.2680434043683275], 'label': 'Dy', 'properties': {}}, {'species': [{'element': 'Dy', 'occu': 1}], 'abc': [0.125, 0.125, 0.125], 'xyz': [1.2680434043683275, 1.2680434043683275, 1.2680434043683275], 'label': 'Dy', 'properties': {}}, {'species': [{'element': 'Ti', 'occu': 1}], 'abc': [0.125, 0.625, 0.625], 'xyz': [6.340217021841638, 3.804130213104983, 3.804130213104983], 'label': 'Ti', 'properties': {}}, {'species': [{'element': 'Ti', 'occu': 1}], 'abc': [0.625, 0.125, 0.625], 'xyz': [3.804130213104983, 6.340217021841638, 3.804130213104983], 'label': 'Ti', 'properties': {}}, {'species': [{'element': 'Ti', 'occu': 1}], 'abc': [0.625, 0.625, 0.125], 'xyz': [3.804130213104983, 3.804130213104983, 6.340217021841638], 'label': 'Ti', 'properties': {}}, {'species': [{'element': 'Ti', 'occu': 1}], 'abc': [0.625, 0.625, 0.625], 'xyz': [6.340217021841638, 6.340217021841638, 6.340217021841638], 'label': 'Ti', 'properties': {}}, {'species': [{'element': 'O', 'occu': 1}], 'abc': [0.7055778297240636, 0.7055778297240636, 0.2944221702759364], 'xyz': [5.07217361747331, 5.07217361747331, 7.157626506000941], 'label': 'O', 'properties': {}}, {'species': [{'element': 'O', 'occu': 1}], 'abc': [0.2944221702759364, 0.2944221702759364, 0.7055778297240636], 'xyz': [5.07217361747331, 5.07217361747331, 2.9867207289456785], 'label': 'O', 'properties': {}}, {'species': [{'element': 'O', 'occu': 1}], 'abc': [0.7055778297240636, 0.2944221702759364, 0.7055778297240636], 'xyz': [5.07217361747331, 7.157626506000941, 5.07217361747331], 'label': 'O', 'properties': {}}, {'species': [{'element': 'O', 'occu': 1}], 'abc': [0.2944221702759364, 0.7055778297240636, 0.2944221702759364], 'xyz': [5.07217361747331, 2.9867207289456785, 5.07217361747331], 'label': 'O', 'properties': {}}, {'species': [{'element': 'O', 'occu': 1}], 'abc': [0.2944221702759364, 0.7055778297240636, 0.7055778297240636], 'xyz': [7.157626506000941, 5.07217361747331, 5.07217361747331], 'label': 'O', 'properties': {}}, {'species': [{'element': 'O', 'occu': 1}], 'abc': [0.7055778297240636, 0.2944221702759364, 0.2944221702759364], 'xyz': [2.9867207289456785, 5.07217361747331, 5.07217361747331], 'label': 'O', 'properties': {}}, {'species': [{'element': 'O', 'occu': 1}], 'abc': [0.5444221702759364, 0.5444221702759364, 0.9555778297240636], 'xyz': [7.608260426209965, 7.608260426209965, 5.522807537682334], 'label': 'O', 'properties': {}}, {'species': [{'element': 'O', 'occu': 1}], 'abc': [0.9555778297240636, 0.9555778297240636, 0.5444221702759364], 'xyz': [7.608260426209965, 7.608260426209965, 9.693713314737597], 'label': 'O', 'properties': {}}, {'species': [{'element': 'O', 'occu': 1}], 'abc': [0.5444221702759364, 0.9555778297240636, 0.5444221702759364], 'xyz': [7.608260426209965, 5.522807537682334, 7.608260426209965], 'label': 'O', 'properties': {}}, {'species': [{'element': 'O', 'occu': 1}], 'abc': [0.9555778297240636, 0.5444221702759364, 0.9555778297240636], 'xyz': [7.608260426209965, 9.693713314737597, 7.608260426209965], 'label': 'O', 'properties': {}}, {'species': [{'element': 'O', 'occu': 1}], 'abc': [0.9555778297240636, 0.5444221702759364, 0.5444221702759364], 'xyz': [5.522807537682334, 7.608260426209965, 7.608260426209965], 'label': 'O', 'properties': {}}, {'species': [{'element': 'O', 'occu': 1}], 'abc': [0.5444221702759364, 0.9555778297240636, 0.9555778297240636], 'xyz': [9.693713314737597, 7.608260426209965, 7.608260426209965], 'label': 'O', 'properties': {}}, {'species': [{'element': 'O', 'occu': 1}], 'abc': [0.25, 0.25, 0.25], 'xyz': [2.536086808736655, 2.536086808736655, 2.536086808736655], 'label': 'O', 'properties': {}}, {'species': [{'element': 'O', 'occu': 1}], 'abc': [0.0, 0.0, 0.0], 'xyz': [0.0, 0.0, 0.0], 'label': 'O', 'properties': {}}]})


if __name__ == "__main__":
    # https://github.com/materialsproject/pymatgen/pull/677
    """
    # need to ask to wiser
    # if structure from file : matrix is ...
    # 0.000000 -5.091292 -5.091292
    # 5.091292 0.000000 -5.091292
    # 5.091292 -5.091292 0.000000
    # and plot relative to axis does not work
    # but if from similar structure : matrix is ...
    # 0.000000 5.072174 5.072174
    # 5.072174 0.000000 5.072174
    # 5.072174 5.072174 0.000000
    # and plot relative to axis does work
    """

    print(structure_pyrochlore_non_magn())
    array = []
    for i in range(22):
        array.append((0, 0, 1))



    # spined_structure = NonCollinearSpinStructure.from_sites_and_orientations(structure_pyrochlore_non_magn(), array, spins_in_sites_or_species="sites")
    # Magmom (spined_structure[2].properties["magmom"])

    from github_non_collinear_spin.SpinIceStructure import ncl_spins_of_typical_structure
    ice_spin = ncl_spins_of_typical_structure.pyrochlore()

    spined_structure = NonCollinearSpinStructure.from_sites_and_orientations(structure_pyrochlore_non_magn(), ice_spin,
                                                                             spins_in_sites_or_species="sites")
    print(spined_structure)

    # vis = VisSpinStruc()
    # vis.set_structure(spined_structure)
    # vis.show()

    # CW = CifWriter(spined_structure, write_magmoms=True)
    # CW.write_file("filename.cif")
