# -*- coding: utf-8 -*-
"""
Description
-----------

Module containing the definition of a DicoVar object, which allow the switch
from different variable denomination in vtk output files, depending on which
solver denomination you use (CharlesX, OpenFOAM, Hybrid, ...)

Contains
--------
dicoVar : _DicoVar object
    The current database for variables denomination
setDicoVar(software) : function
    Set the current database to a given software denomination
printDatabase() : function
    Print the current database denomination

Examples
--------
When implementing a vtk file manipulation function

>>> # Import dicovar module database
>>> from aerokit.common.dicovar import dicoVar as dv
>>> # test if the field U_AVG is present
>>> if data_outVTK.GetPointData().HasArray(dv('u_avg')) != 1:
>>>     raise ValueError("Error : field U_AVG not present")
>>> # test if the field RHO_AVG is present
>>> if data_outVTK.GetPointData().HasArray(dv['rho_avg']) != 1:
>>>     raise ValueError("Error : field RHO_AVG not present")
>>> # test if the field TauWallAvg is present
>>> if wall_outVTK.GetPointData().HasArray(dv.get('tauw_avg')) != 1:
        raise ValueError("Error : field TauWallAvg not present")

When implementing a script associated to a given software

>>> # Import dicovar module functions
>>> from aerokit.common.dicovar import setDicoVar, printDicoVar
>>> # Setting CharlesX dicovar
>>> setDicoVar('cx')
>>> printDicoVar()

Notes
-----
- The implementation allow to use several names to set one software (ex : 'cx',
  'CharlesX', ... for CharlesX denomination)
- Several ways to get a variables : dv[...], dv(...), dv.get(...)
- variable name handles capital letters : dv['U_avg'] <=> dv['u_avg']
"""
import copy as _copy

# -----------------------------------------------------------------
# Part to modify in order to implement other variable denominations
# ---------------

# CharlesX variables dico
_cxNames = ['IC3', 'ic3', 'CharlesX', 'cx', 'charlesx']
_cxDico = {
    'P':   'P',
    'p':   'P',
    'Ps':  'P',
    'rho': 'RHO',
    'T':   'T',
    'V':  'U',
    'Vx': 'U_X',
    'Vy': 'U_Y',
    'Vz': 'U_Z',
    'u_avg': 'U_AVG',
    'rho_avg': 'RHO_AVG',
    'tauw_avg': 'TauWallAvg',
    'mulam_avg': 'MU_LAM_AVG',
    't_avg': 'T_AVG',
    'u_rms': 'U_RMS',
    'u_rey': 'U_REY'}

# OpenFOAM variables dico (just examples, to modify ...)
_ofNames = ['OpenFOAM', 'of', 'openfoam', 'ofoam']
_ofDico = {
    'u_avg': 'Uavg',
    'rho_avg': 'rhoAVG',
    'tauw_avg': '???',
    'mulam_avg': '???',
    't_avg': '???',
    'u_rms': '???',
    'u_rey': '???'}

# Database with all software denominations
_dataBase = [[_cxNames, _cxDico],
             [_ofNames, _ofDico]]

_default = 'IC3'


# -----------------------------------------------------------------
# DO NOT MODIFY THIS PART ...
# -----------------------

class _DicoVar(object):
    """
    Class representing a database for variables denomination of CFD
    softwares
    """

    def __init__(self, software):

        # Internal variables definitions
        self._dv = {}
        self._softName = ''

        # Setting database
        softOK = False
        for soft in _dataBase:
            if software.lower() in soft[0]:
                self._dv = _copy.deepcopy(soft[1])
                self._softName = soft[0][0]
                softOK = True
                print('---- HADES::DicoVar ----')
                print('--> Using {} variables denomination'.format(self._softName))
        if not softOK:
            raise ValueError('Not variable denomination database for {}'
                             .format(software))

    # Main class methods
    def keys(self):
        """
        return list of keys in dictionnary
        """
        return self._dv.keys()            

    def get(self, var):
        """
        Get the variable name with its denomination in the selected software database
        Parameters:
            var : String            The DaePy denomination of the variable
        Returns:
            varSoftware : string            The selected software denomination of the variable
        """
        if var in self._dv:
            return self._dv[var]
        else:
            raise ValueError('{} variable not in {} database dictionnary'.format(var, self._softName))

    def addVar(self, var, sVar):
        """
        Add a variable name with its denomination in the selected software database
        Parameters :
            var : String        The DaePy denomination of the variable (not in database)
            sVar : string       The software denomination of the variable
        """
        if var not in self._dv:
            self._dv[var] = sVar
        else:
            raise ValueError('{} variable already in {}'.format(var, self._softName) + 'database dictionnary.\n' +'Use setVar function instead')
        
    def setVar(self, var, sVar):
        """
        Set a variable name with a new denomination in the selected software database
        Parameters
        ----------
        var : String            The DaePy denomination of the variable (already in database)
        sVar : String           The new software denomination of the variable
        """
        if var in self._dv:
            self._dv[var] = sVar
            print('---- HADES::DicoVar ----')
            print('--> Changing denomination of {} into {}'.format(var, sVar))
        else:
            raise ValueError('{} variable not in {} database dictionnary'.format(var, self._softName))

    def printDatabase(self):
        """Print the variables denomination for the selected software"""
        print('---- HADES::DicoVar ----')
        print('--> Database dictionnary for {}'.format(self._softName))
        for elt in self._dv:
            print('      {} => {}'.format(elt, self._dv[elt]))

    # Definition of other ways to get a variable
    def __call__(self, var):
        return self.get(var)

    def __getitem__(self, var):
        return self.get(var)




# The main _DicoVar variable
dicoVar = _DicoVar(_default)


def setDicoVar(software):
    """
    Associate the denomination database to a given software

    Parameters
    ----------
    software : String
        The selected software (CharlesX, OpenFOAM, Hybrid, ...)
    """
    dicoVar.__init__(software)


def printDicoVar():
    """Print the variables denominations for the selected software"""
    dicoVar.printDatabase()

# -----------------------------------------------------------------

