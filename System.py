'''
@file
@date 2014
@author Karl R. Leikanger
'''
class System():
    '''
    @brief
    @var
    @date 2014
    @author Karl R. Leikanger
    '''
    
    class atom():
        '''
        @brief
        @var basisset Basisset for atom.
        @var position Position of atom in the reference cell. 
        @var ghost Is it a ghost?
        @date 2014
        @author Karl R. Leikanger
        '''
        basisset = ''
        position = ''
        ghost = False
    
    class lattice():
        '''
        @brief
        @var avec Lattice vectors.
        @var nlattice Integer array size of supercell in terms of avec.
        @date 2014
        @author Karl R. Leikanger
        '''
        avec = []
        nlattice = []

    def __init__()
       '''
       @brief
       @date 2014
       @author Karl R. Leikanger
       '''

    def read_initfile(self):
        '''
        @brief
        @param
        @date 2014
        @author Karl R. Leikanger
        '''

    def write_cp2k_inputfile()
        '''
        @brief
        @param
        @date 2014
        @author Karl R. Leikanger
        '''



