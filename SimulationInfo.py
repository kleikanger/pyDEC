'''
@file Contains a class with info about the current simulations.
@brief
@param
@date 2014
@author Karl R. Leikanger

Parameters like inputfiles, the basis, latticevectors etc are stored in this
class.
'''

# import basis
import parameters
from IO import ReadFile
# import inputparam
import sys
import System


class SimulationInfo():
    '''
    @brief
    @var
    @date 2014
    @author Karl R. Leikanger
    '''
    basis = ''
    system = ''
    # inputparam = ''

    def __init__(self, basisname, basisfile, basisformat, sysfile, inputfile):
        '''
        @brief
        @param
        @date 2014
        @author Karl R. Leikanger
        '''
        # self.basis = Basis(basisname)
        # self.basis.read_dalton_basisfile(basisfile) #basisfile, basis_format

        # self.inputpatam = InputParam()
        # self.InputParam.read(inputfile)
        # cp2k default settings file
        # dalton default settings file

        self.system = System()
        self.system.read(sysfile)

