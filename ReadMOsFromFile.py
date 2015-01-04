# see CP2K cp_fm_types.F
import IO

try:
    import numpy as np
except:
    print('\nError loading nympy library in python.')
    print('Make sure that the library is available.\n')
    raise

'''
@brief Read restart file in CP2K format and convert to LSDALTON format.
@param
@date 2014
@author Karl R. Leikanger.
'''


class MOinfo():
    '''
    @brief
    @param
    @date 2014
    @author Karl R. Leikanger.
    '''

    filename = ''
    int_type = ''
    double_type = ''
    rt_mos = ''         # TODO should be input
    i_nmo = 9999999999  # =i_ao # TODO should be input ?
    # set from basis_set : sum of max(ao.tag) for all elem/ basis sets.

    '''
    basis_name = ''
    nspin = 0
    nshell = []
    l = []
    nshell_info = []
    nso_info = []

    ...
    '''

    def __init__(self, filename, rt_mos=False,
                 int_type=np.int32, double_type=np.float64):
        '''
        @brief The constructor also reads the CP2K datafile (filename)
        @brief and loads the data from the CP2K simulation.
        @param
        @date 2014
        @author Karl R. Leikanger.
        '''
        self.rt_mos = rt_mos
        self.filename = filename
        self.int_type = int_type
        self.double_type = double_type

    def read_restart_file_cp2k_format(self, system, basis):  # convert_mo_format
        '''
        @brief Load CP2K input file.
        @param system xxx
        @param basis xxx
        @date 2014
        @author Karl R. Leikanger.

        See CP2K:qs_mo_io:read_mos_restart_low

        For this routine to work, must the dimension <number of ao's>
        be the same in both simulations ??
        '''
        read_record = IO.ReadFortranBinaryFile(
            self.int_type, self.double_type, self.filename
        ).read_record

        # Read general info.
        [i_natom, i_nspin, i_nao, i_nsetmax, i_nshellmax]\
            = read_record('INT', 5)
        print([i_natom, i_nspin, i_nao, i_nsetmax, i_nshellmax])
        nsetinfo = read_record('INT', i_natom)
        nshell_info = read_record('INT', i_nsetmax*i_natom)
        nshell_info.reshape(i_nsetmax, i_natom)
        nso_info = read_record('INT', i_nshellmax*i_nsetmax*i_natom)
        nso_info.reshape(i_nshellmax, i_nsetmax, i_natom)

        # read MO's
        i_nmo = self.i_nmo
        for ispin in range(1, i_nspin+1):

            # IF (para_env%ionode.AND.(nmo > 0)) THEN
            [i_nmo_read, i_homo, i_lfomo, i_nelectron] =\
                read_record('INT', 4)

            # TODO This might not be neccessary. Check LSDALTON/CP2K source.
            # But note that in LSDALTON it is assured that the allocated matrix
            # have the correct ncol and nrow (nmo and nao).
            i_nmo = min(i_nmo, i_nmo_read)
            if (i_nmo_read > i_nmo):
                print('The number of MOs on the restart file is smaller than\
                      the number of the allocated MOs in LSDALTON. The MO set\
                      will be padded with zeros!')
            elif (i_nmo_read < i_nmo):
                print(' The number of MOs on the restart file is greater\
                      than the number of the allocated MOs in LSDALTON.\
                      The read MO set will be truncated!')
            #
            # Various tests to check that nmo > i_homo etc.
            # ...
            #

            # Read and pad if i_nmo_read>i_nmo
            # XXX sorted? If so: last eigs have largest e and will be truncated
            tmp = read_record('DOUBLE', i_nmo_read*2)
            dd_eig = tmp[0:i_nmo]
            dd_occ = tmp[i_nmo_read:i_nmo_read+i_nmo]
            # if self.i_nmo > i_nmo_read:
            #    dd_eig += [0] * (self.i_nmo - i_nmo_read)
            #    dd_occ += [0] * (self.i_nmo - i_nmo_read)
            # XXX Use when self.i_nmo is set to the correct value.
            # XXX these variables (dd_eig, dd_occ) are not used.

            permutations =\
                basis.get_mo_transformation('CP2K', 'DALTON', system.atoms)

            # open output filestream
            ostream = IO.WriteFortranBinaryFile(
                self.int_type, self.double_type, 'orbitals_in.u'
            )

            nrow = i_nao
            ncol = i_nmo

            ostream.write_record('INT', [i_nao, self.i_nmo])
            # Buffer: start of MO - matrix. Enables us to write
            # MO's one line at a time.
            ostream.write_record('INT', nrow*ncol, buf=False)

            if self.rt_mos:  # unres or res??
                print('Error ReadMOsFromFile: Support for unrestricted MOs\
                      not implemented.')
                raise SystemError

                for imat in range(2*ispin-1, 2*ispin+1):
                    for i in range(1, i_nmo+1):
                        dd_vecbuff = read_record('DOUBLE', i_nao)
                        '''
                        see in qs_mo_io.F what to do with the data
                        '''
            else:
                for i in range(i_nmo):

                    # NB: this is a numpy matrix.
                    dd_vecbuff = read_record('DOUBLE', i_nao)
                    print(dd_vecbuff)
                    #dd_vecbuff = [dd_vecbuff[j] for j in permutations]
                    dd_vecbuff = dd_vecbuff[permutations]
                    print(dd_vecbuff)

                    ostream.write_record('DOUBLE', dd_vecbuff, buf=False)

            '''
            # XXX only necc. if i_nspin > 1
            # read the rest of the MOs (if there are any)
            if i_nmo < i_nmo_read:
                n = i_nmo_read - i_nmo
                for i in range(n)
                    read_record('DOUBLE', i_nao)
            # or add buffer of zeroes
            if self.i_nmo > i_nmo_read:
                buf = np.zeros(i_ao)
                n = self.i_nmo - i_nmo_read
                for i in range(n):
                    ostream.write_record('DOUBLE', buf, buf=False)

            '''
            # XXX Use above when self.i_nmo is set to the correct value.

            # Buffer: end of MO - matrix.
            ostream.write_record('INT', nrow*ncol, buf=False)





from System import System
system = System()
system.read_system_info('N2.mol', 'LSDALTON')

# o = MOinfo('He_bulk2-RESTART.wfn')
# # convert_restart_file ,,
# o.read_restart_file_cp2k_format(system, basis)

from basis import Basis
basis = Basis()
basis.read_basis_set('STO-3G', 'STO-3G', 'DALTON')
basis.print_basis_sets('cp2k_STO-3G__3', 'CP2K')

o = MOinfo('Ne_bulk2-RESTART.wfn')
o.read_restart_file_cp2k_format(system, basis)
print('ao order DALTON')
print(basis.get_mo_transformation('CP2K', 'DALTON', system.atoms))
