from github_non_collinear_spin.NonColSpinStructure import NonCollinearSpinStructure
from typing import Dict, List, Tuple, Optional, Union, Iterator, Set, Sequence, Iterable
from pymatgen.core.lattice import Lattice, get_points_in_spheres
from pymatgen.core.periodic_table import Element, Specie, get_el_sp, DummySpecie
from pymatgen.core.composition import Composition
from pymatgen.core.structure import Structure
from pymatgen.symmetry.analyzer import SpacegroupAnalyzer
from pymatgen.vis.structure_vtk import StructureVis
import vtk
from pymatgen.electronic_structure.core import Magmom
from pymatgen.core.periodic_table import Element
from pymatgen.core.periodic_table import Specie
import numpy as np
from pymatgen.core.sites import PeriodicSite
import os
from monty.serialization import loadfn, dumpfn
import pymatgen.analysis.magnetism.analyzer
import numpy as np
from pymatgen.io.cif import CifWriter
from collections import OrderedDict, deque
from ast import literal_eval

MODULE_DIR = os.path.dirname(os.path.abspath(pymatgen.analysis.magnetism.analyzer.__file__))
DEFAULT_MAGMOMS = loadfn(os.path.join(MODULE_DIR, "default_magmoms.yaml"))



def structure_pyrochlore_non_magn():
    return Structure.from_dict({'@module': 'pymatgen.core.structure', '@class': 'Structure', 'charge': None, 'lattice': {'matrix': [[0.0, 5.07217361747331, 5.07217361747331], [5.07217361747331, 0.0, 5.07217361747331], [5.07217361747331, 5.07217361747331, 0.0]], 'a': 7.173136720541758, 'b': 7.173136720541758, 'c': 7.173136720541758, 'alpha': 60.00000000000001, 'beta': 60.00000000000001, 'gamma': 60.00000000000001, 'volume': 260.9830654620022}, 'sites': [{'species': [{'element': 'Dy', 'occu': 1}], 'abc': [0.625, 0.125, 0.125], 'xyz': [1.2680434043683275, 3.804130213104983, 3.804130213104983], 'label': 'Dy', 'properties': {}}, {'species': [{'element': 'Dy', 'occu': 1}], 'abc': [0.125, 0.625, 0.125], 'xyz': [3.804130213104983, 1.2680434043683275, 3.804130213104983], 'label': 'Dy', 'properties': {}}, {'species': [{'element': 'Dy', 'occu': 1}], 'abc': [0.125, 0.125, 0.625], 'xyz': [3.804130213104983, 3.804130213104983, 1.2680434043683275], 'label': 'Dy', 'properties': {}}, {'species': [{'element': 'Dy', 'occu': 1}], 'abc': [0.125, 0.125, 0.125], 'xyz': [1.2680434043683275, 1.2680434043683275, 1.2680434043683275], 'label': 'Dy', 'properties': {}}, {'species': [{'element': 'Ti', 'occu': 1}], 'abc': [0.125, 0.625, 0.625], 'xyz': [6.340217021841638, 3.804130213104983, 3.804130213104983], 'label': 'Ti', 'properties': {}}, {'species': [{'element': 'Ti', 'occu': 1}], 'abc': [0.625, 0.125, 0.625], 'xyz': [3.804130213104983, 6.340217021841638, 3.804130213104983], 'label': 'Ti', 'properties': {}}, {'species': [{'element': 'Ti', 'occu': 1}], 'abc': [0.625, 0.625, 0.125], 'xyz': [3.804130213104983, 3.804130213104983, 6.340217021841638], 'label': 'Ti', 'properties': {}}, {'species': [{'element': 'Ti', 'occu': 1}], 'abc': [0.625, 0.625, 0.625], 'xyz': [6.340217021841638, 6.340217021841638, 6.340217021841638], 'label': 'Ti', 'properties': {}}, {'species': [{'element': 'O', 'occu': 1}], 'abc': [0.7055778297240636, 0.7055778297240636, 0.2944221702759364], 'xyz': [5.07217361747331, 5.07217361747331, 7.157626506000941], 'label': 'O', 'properties': {}}, {'species': [{'element': 'O', 'occu': 1}], 'abc': [0.2944221702759364, 0.2944221702759364, 0.7055778297240636], 'xyz': [5.07217361747331, 5.07217361747331, 2.9867207289456785], 'label': 'O', 'properties': {}}, {'species': [{'element': 'O', 'occu': 1}], 'abc': [0.7055778297240636, 0.2944221702759364, 0.7055778297240636], 'xyz': [5.07217361747331, 7.157626506000941, 5.07217361747331], 'label': 'O', 'properties': {}}, {'species': [{'element': 'O', 'occu': 1}], 'abc': [0.2944221702759364, 0.7055778297240636, 0.2944221702759364], 'xyz': [5.07217361747331, 2.9867207289456785, 5.07217361747331], 'label': 'O', 'properties': {}}, {'species': [{'element': 'O', 'occu': 1}], 'abc': [0.2944221702759364, 0.7055778297240636, 0.7055778297240636], 'xyz': [7.157626506000941, 5.07217361747331, 5.07217361747331], 'label': 'O', 'properties': {}}, {'species': [{'element': 'O', 'occu': 1}], 'abc': [0.7055778297240636, 0.2944221702759364, 0.2944221702759364], 'xyz': [2.9867207289456785, 5.07217361747331, 5.07217361747331], 'label': 'O', 'properties': {}}, {'species': [{'element': 'O', 'occu': 1}], 'abc': [0.5444221702759364, 0.5444221702759364, 0.9555778297240636], 'xyz': [7.608260426209965, 7.608260426209965, 5.522807537682334], 'label': 'O', 'properties': {}}, {'species': [{'element': 'O', 'occu': 1}], 'abc': [0.9555778297240636, 0.9555778297240636, 0.5444221702759364], 'xyz': [7.608260426209965, 7.608260426209965, 9.693713314737597], 'label': 'O', 'properties': {}}, {'species': [{'element': 'O', 'occu': 1}], 'abc': [0.5444221702759364, 0.9555778297240636, 0.5444221702759364], 'xyz': [7.608260426209965, 5.522807537682334, 7.608260426209965], 'label': 'O', 'properties': {}}, {'species': [{'element': 'O', 'occu': 1}], 'abc': [0.9555778297240636, 0.5444221702759364, 0.9555778297240636], 'xyz': [7.608260426209965, 9.693713314737597, 7.608260426209965], 'label': 'O', 'properties': {}}, {'species': [{'element': 'O', 'occu': 1}], 'abc': [0.9555778297240636, 0.5444221702759364, 0.5444221702759364], 'xyz': [5.522807537682334, 7.608260426209965, 7.608260426209965], 'label': 'O', 'properties': {}}, {'species': [{'element': 'O', 'occu': 1}], 'abc': [0.5444221702759364, 0.9555778297240636, 0.9555778297240636], 'xyz': [9.693713314737597, 7.608260426209965, 7.608260426209965], 'label': 'O', 'properties': {}}, {'species': [{'element': 'O', 'occu': 1}], 'abc': [0.25, 0.25, 0.25], 'xyz': [2.536086808736655, 2.536086808736655, 2.536086808736655], 'label': 'O', 'properties': {}}, {'species': [{'element': 'O', 'occu': 1}], 'abc': [0.0, 0.0, 0.0], 'xyz': [0.0, 0.0, 0.0], 'label': 'O', 'properties': {}}]})



class ncl_spins_of_typical_structure():

    @staticmethod
    def pyrochlore(site_1='in', site_2='in', site_3='in', site_4='in', user_magmom=False):

        dict_spin_ice_prim_cell = {'[0.625, 0.125, 0.125]': {'in': (1, -1, -1), 'out': [-1, 1, 1]},
                                   '[0.125, 0.625, 0.125]': {'in': [-1, 1, -1], 'out': [1, -1, 1]},
                                   '[0.125, 0.125, 0.625]':{'out': [1, 1, -1],'in': [-1, -1, 1]},
                                   '[0.125, 0.125, 0.125]':{'in': [1, 1, 1],'out': [-1, -1, -1]}}

        species_A = ['[0.625, 0.125, 0.125]', '[0.125, 0.625, 0.125]', '[0.125, 0.125, 0.625]', '[0.125, 0.125, 0.125]']
        species_B = ['[0.125, 0.625, 0.625]', '[0.625, 0.125, 0.625]', '[0.625, 0.625, 0.125]', '[0.625, 0.625, 0.625]']
        species_O = ['O' * 14]

        array_magn = []
        array_magn.append(dict_spin_ice_prim_cell['[0.625, 0.125, 0.125]'][site_1])
        array_magn.append(dict_spin_ice_prim_cell['[0.125, 0.625, 0.125]'][site_2])
        array_magn.append(dict_spin_ice_prim_cell['[0.125, 0.125, 0.625]'][site_3])
        array_magn.append(dict_spin_ice_prim_cell['[0.125, 0.125, 0.125]'][site_4])
        for i in range(18):
            array_magn.append([0, 0, 1])
        # d_magn = {}
        # d_magn['[0.625, 0.125, 0.125]'] = dict_spin_ice_prim_cell['[0.625, 0.125, 0.125]'][site_1]
        # d_magn['[0.125, 0.625, 0.125]'] = dict_spin_ice_prim_cell['[0.125, 0.625, 0.125]'][site_2]
        # d_magn['[0.125, 0.125, 0.625]'] = dict_spin_ice_prim_cell['[0.125, 0.125, 0.625]'][site_3]
        # d_magn['[0.125, 0.125, 0.125]'] = dict_spin_ice_prim_cell['[0.125, 0.125, 0.125]'][site_4]
        return array_magn



print(ncl_spins_of_typical_structure.pyrochlore())
