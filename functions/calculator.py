"""
The module contains a function 'get_calc' that creates a calculator SIESTA with given parameters.
"""

from ase.calculators.siesta import Siesta
from ase.units import Ry

def get_calc(label):
    """
    Returns the calculator object for the given label.

    Args:
        label (str): The label to identify output files.

    Returns:
        ase.Calculator: The calculator object with given parameters.
    """
    tmp_calc = Siesta(label=f'{label}',
                  xc='PBE',
                  basis_set='SZP',
                  kpts=[8, 8, 1],
                  mesh_cutoff=200 * Ry,
                  spin='collinear',
                  fdf_arguments={'MaxSCFIterations':    500,
                                'PAO.BasisSize':           'SZP',
                                'MD.NumCGsteps':        250,
                                'MD.TypeOfRun':         'CG',
                                'MD.VariableCell':      'F',
                                'MD.MaxCGDispl':        '0.2000000000  Bohr',
                                'MD.MaxForceTol':       '0.1 eV/Ang',
                                'SolutionMethod':     'Diagon',
                                'DM.MixingWeight':   0.1000000000,
                                'DM.NumberPulay':     '6',
                                'DM.NumberBroyden':        0,
                                'DM.OccupancyTolerance':   0.1000000000E-11,
                                'DM.NumberKick':           0,
                                'DM.KickMixingWeight':     0.5000000000,
                                'DM.Tolerance':       1.0E-3,
                                'DM.UseSaveDM':       '.true.',
                                'DM.UseSaveXV':       '.true.',
                                'WriteMullikenPop':        1 ,
                                'WriteDenchar':            '.true.',
                                'WriteKpoints':            '.true.',
                                'WriteForces':             '.true.',
                                'WriteDM':                 '.true.',
                                'WriteXML':                '.true.',
                                'WriteEigenvalues':        '.false.',
                                'WriteCoorStep':           '.true.',
                                'WriteMDhistory':          '.true.',
                                'WriteMDXmol':             '.true.',
                                'WriteCoorXmol':           '.true.',
                                 },
                 )
    return tmp_calc