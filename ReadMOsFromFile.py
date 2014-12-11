#see CP2K cp_fm_types.F
import sys
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
    rt_mos = ''        #TODO should be input
    i_nmo = 9999999999 #TODO should be input ?
    
    '''
    basis_name = ''
    nspin = 0
    nshell = []
    l = []
    nshell_info = []
    nso_info = []

    ...
    '''

    def __init__(self, filename, \
            rt_mos = False, \
            int_type = np.int32, \
            double_type = np.float64):
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

        self.read_restart_file_cp2k_format()


    def read_record(self, f, dtype, nread):
        '''
        @brief Read a single record from a fortran output file.
        @param f Inputstream.
        @param datatype Numpy datatype.
        @param nread Number of datatype elements to read. 
        @return retvec Numpy array with read data.
        @date 2014
        @author Karl R. Leikanger.

        Fortran IO is record based, not stream based:
        For unformated IO, Fortran compilers typically write the length of 
        the record at the beginning and end of the record. Most but not all 
        compilers use four bytes. 

        '''

        buf_type = self.int_type

        c1 = np.fromfile(f, dtype=buf_type, count=1)[0]
        retvec = np.fromfile(f, dtype=dtype, count=nread)
        c2 = np.fromfile(f, dtype=buf_type, count=1)[0]

        if (c1 != c2):
            #TODO Detailed error message
            print('\nError in input file. Check output format, datatypes etc.\n')
            exit(-1)

        return retvec

    def read_restart_file_cp2k_format(self): 
        '''
        @brief Load CP2K input file.
        @date 2014
        @author Karl R. Leikanger.

        '''

        int_type = self.int_type
        double_type = self.double_type
        int_nbytes = int_type(0).nbytes
        double_nbytes = double_type(0).nbytes

        print('Reading binary file %s assuming %iB integer and %iB float.' \
                %(self.filename, int_nbytes, double_nbytes))

        try:
            f = open(self.filename, "rb")
        except IOError as e:
            print("Error while trying to open %s", self.filename)
            print("IO error({0}): {1}".format(e.errno, e.strerror))
        except:
            print('unexcepted error', sys.exc_info()[0])
            raise

        # Read the file. See CP2K:qs_mo_io:read_mos_restart_low

        [i_natom, i_nspin, i_nao, i_nsetmax, i_nshellmax] =\
                self.read_record(f, int_type, 5)

        nsetinfo = self.read_record(f, int_type, i_natom)

        nshell_info = self.read_record(f, int_type, i_nsetmax*i_natom)
        nshell_info.reshape(i_nsetmax, i_natom)

        nso_info = self.read_record(f, int_type, i_nshellmax*i_nsetmax*i_natom)
        nso_info.reshape(i_nshellmax, i_nsetmax, i_natom)

        i_nmo = self.i_nmo

        for ispin in range(1, i_nspin+1):
            #IF (para_env%ionode.AND.(nmo > 0)) THEN
            [i_nmo_read, i_homo, i_lfomo, i_nelectron] =\
                    self.read_record(f, int_type, 4)

            #TODO assert (nmo_read >= nmo)
            i_nmo = min(i_nmo, i_nmo_read) #??
            
            dd_eig = np.zeros(i_nmo, dtype=double_type)
            dd_occ = np.zeros(i_nmo, dtype=double_type)
            tmp =  self.read_record(f, double_type, i_nmo_read*2)
            dd_eig[0:i_nmo] = tmp[0:i_nmo]  
            dd_occ[0:i_nmo] = tmp[i_nmo_read:i_nmo_read+i_nmo]

            if self.rt_mos: #unres or res??
                for imat in range(2*ispin-1, 2*ispin+1):	
                    for i in range(1, i_nmo+1):

                        dd_vecbuff = self.read_record(f, double_type, i_nao)
                        '''
                        see in qs_mo_io.F what to do with the data
                        '''
            else:
                for i in range(1, i_nmo+1):

                    dd_vecbuff = self.read_record(f, double_type, i_nao)


                    #
                    # Input info
                    #
                    '''
                    for iatom in range(1, natom+1):
                        get_atomic_kind
                        get_qs_kind (??)
                        if associated(orb_basis_set) then
                            get_gto_basis
                            minbas=false
                        elseif associated(dftb_parameter) then
                            get_dftb_atom_param
                            minbas=false
                        else
                            unknown basis set?

                         IF (para_env%ionode) THEN

                    '''
                    #
                    # Reorganize vecbuf
                    #
                    '''

                  use_this = .TRUE.
                  iset_read = 1
                  DO iset=1,nset 
                  '''
                  # Input from basis file. 
                  # the number of sets in basis, each contains 
                  # nexp * (lshell(lmin) + ... + nshell(lmax)) ao's
                  '''
                     ishell_read = 1
                     IF(minbas) THEN
                        nnshell = lmax+1
                     ELSE
                        nnshell = nshell(iset)
                        '''
                        # Input from basis file. 
                        # nshell input
                        '''
                     END IF
                     DO ishell=1,nnshell
                        IF(minbas) THEN
                           lshell = ishell-1
                        ELSE
                           lshell = l(ishell,iset)
                        END IF
                        IF (iset_read > nset_info(iatom)) use_this = .FALSE.
                        IF (use_this) THEN ! avoids out of bound access of the lower line if false
                           IF (nso(lshell) == nso_info(ishell_read,iset_read,iatom)) THEN
                              offset_read=offset_info(ishell_read,iset_read,iatom)
                              ishell_read=ishell_read+1
                              IF (ishell_read > nshell_info(iset,iatom)) THEN
                                 ishell_read = 1
                                 iset_read = iset_read+1
                              END IF
                           ELSE
                              use_this = .FALSE.
                           END IF
                        END IF
                        DO iso=1,nso(lshell)
                           IF (use_this) THEN
                              IF (offset_read-1+iso.LT.1 .OR. offset_read-1+iso.GT.nao_read) THEN
                                 vecbuffer(1,irow)=0.0_dp
                              ELSE
                                 vecbuffer(1,irow)=vecbuffer_read(1,offset_read-1+iso)
                              END IF
                           ELSE
                              vecbuffer(1,irow) = 0.0_dp
                           END IF
                           irow = irow + 1
                        END DO
                        use_this = .TRUE.
                        '''
                        #then set MO's
                        '''
                CALL mp_bcast(vecbuffer,source,group)
                CALL cp_fm_set_submatrix(mos(ispin)%mo_set%mo_coeff,&
                 vecbuffer,1,i,nao,1,transpose=.TRUE.,error=error)

            '''
            #Read the rest of the nmo's if there are any...
            '''
            IF(nmo>0) THEN
               DO i=nmo+1,nmo_read
                    READ (rst_unit) vecbuffer_read

                    '''

        f.close()
''' 

         #      IF (para_env%ionode) THEN
         #         READ (rst_unit) vecbuffer
         #      ELSE
         #         vecbuffer(1,:) = 0.0_dp
         #      END IF
         #      CALL mp_bcast(vecbuffer,source,group)
         #      CALL cp_fm_set_submatrix(rt_mos(imat)%matrix,&
         #           vecbuffer,1,i,nao,1,transpose=.TRUE.,error=error)
         #   END DO
         #END DO
         #ELSE
			#  DO i=1,nmo
         #   IF (para_env%ionode) THEN
         #      READ (rst_unit) vecbuffer_read
         #      ! now, try to assign the read to the real vector
         #      ! in case the basis set changed this involves some guessing
         #      irow=1
         #      DO iatom=1,natom
			#
			#
			#! Skip extra MOs if there any
'''		


#	def print_restart_file_dalton_format(self):

#		! *****************************************************************************




#__debug__ set to?
#assert a<b, 'errorstring'
moinfo = MOinfo('slettmeg.wfn')




'''

Below a function to read mo coeff from unit ires

'''

#!> \brief ...
#!> \param mo_set ...
#!> \param ires ...
#!> \param error ...
#! *****************************************************************************
#  SUBROUTINE read_mo_set_basic (mo_set, ires, error)
#
#    TYPE(mo_set_type), POINTER               :: mo_set
#    INTEGER, INTENT(in)                      :: ires
#    TYPE(cp_error_type), INTENT(inout)       :: error
#
#    CHARACTER(LEN=*), PARAMETER :: routineN = 'read_mo_set_basic', &
#      routineP = moduleN//':'//routineN
#
#    INTEGER                                  :: mao, mmo, nao, nmo
#    LOGICAL                                  :: failure
#
#    mao = mo_set%nao
#    mmo = mo_set%nmo
#
#    IF ((ires>0).AND.(nmo > 0)) THEN
#       READ (ires) nao,nmo,mo_set%homo,mo_set%lfomo,mo_set%nelectron
#       CALL cp_assert(mao==nao,cp_fatal_level,cp_assertion_failed,routineP,&
#                         "Different size of basis set",error,failure)
#       CALL cp_assert(mmo==nmo,cp_fatal_level,cp_assertion_failed,routineP,&
#                         "Different number of MOs",error,failure)
#       READ (ires) mo_set%eigenvalues(1:nmo),mo_set%occupation_numbers(1:nmo)
#    END IF
#    CALL cp_fm_read_unformatted(mo_set%mo_coeff,ires,error)
#
#  END SUBROUTINE read_mo_set_basic

'''
Read unformatted is 
'''

#    CALL cp_fm_read_unformatted(mo_set%mo_coeff,ires,error)
#to
#SUBROUTINE cp_fm_read_unformatted(fm,unit,error)
#    TYPE(cp_fm_type), POINTER                :: fm
#...
#...
#
#   DO j=1,ncol_global
#      READ (unit) fm%local_data(:,j)
#   END DO
'''
The struct looks like this
'''

#  TYPE cp_fm_type
#!    PRIVATE
#     CHARACTER(LEN=60) :: name
#     INTEGER :: id_nr, ref_count, print_count
#     LOGICAL :: use_sp
#     TYPE(cp_fm_struct_type), POINTER :: matrix_struct
#     REAL(KIND = dp), DIMENSION(:,:), POINTER :: local_data
#     REAL(KIND = sp), DIMENSION(:,:), POINTER :: local_data_sp
#  END TYPE cp_fm_type


'''
        To complicate the situation: some compilers such as Gfortran and Intel
        Fortran support records larger than 2 GB despite having 4 byte record 
        markers, by using subrecords.

        This is the reason why it may be smarter to convert the files using
        a fortran routine that is compiled with the same compiler as CP2K.

        Here, we simply use 4B buffers. But obvuously this is not a optimal 
        solution.
'''

''' 
Basis set format: CP2K
======================
Element symbol  Name of the basis set  Alias names
nset (repeat the following block of lines nset times)
n lmin lmax nexp nshell(lmin) nshell(lmin+1) ... nshell(lmax-1) nshell(lmax)
a(1)      c(1,l,1)      c(1,l,2) ...      c(1,l,nshell(l)-1)      c(1,l,nshell(l)), l=lmin,lmax
a(2)      c(2,l,1)      c(2,l,2) ...      c(2,l,nshell(l)-1)      c(2,l,nshell(l)), l=lmin,lmax
 .         .             .                 .                       .
 .         .             .                 .                       .
 .         .             .                 .                       .
a(nexp-1) c(nexp-1,l,1) c(nexp-1,l,2) ... c(nexp-1,l,nshell(l)-1) c(nexp-1,l,nshell(l)), l=lmin,lmax
a(nexp)   c(nexp,l,1)   c(nexp,l,2)   ... c(nexp,l,nshell(l)-1)   c(nexp,l,nshell(l)), l=lmin,lmax


nset     : Number of exponent sets
n        : Principle quantum number (only for orbital label printing)
lmax     : Maximum angular momentum quantum number l
lmin     : Minimum angular momentum quantum number l
nshell(l): Number of shells for angular momentum quantum number l
a        : Exponent
c        : Contraction coefficient
'''
