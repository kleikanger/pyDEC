
'''
@file Contains a class that reads and writes basisfiles in both CP2K and DALTON format.
@brief Read basis file, store basis info and generate basis input to CP2K or DALTON
@date 2014
@author Karl R. Leikanger.

import sys
import parameters as param
from IO import ReadFile

class Basis():
    '''
    @brief
    @var chem_elems list of class ChemElemBasis objects. Basis for given atomic nr. 
    @var chem_elems_indx dictionary {atomic_nr : index in chem_elems }
    @var basis_name Name of basis.
    @date 2014
    @author Karl R. Leikanger.
    '''

    chem_elems = [] 
    chem_elems_indx = '' 
    basis_name = ''
    
    class ChemElemBasis():
        '''
        @brief Container class. Stores basis info for a given element.
        @var atomic_nr Atomic number.
        @var ao List of AtomicOrbital objects.
        @date 2014
        @author Karl R. Leikanger.
        '''
        atomic_nr = '' 
        ao = []
    
    class AtomicOrbital():
        '''
        @brief Container class. Stores basis info for a given element.
        @var exp Exponents
        @var coef Contraction coefficients
        @var l Quantum number L.
        @date 2014
        @author Karl R. Leikanger.
        '''
        e = [] 
        c = []
        l = ''
        
    def __init__(self, basis_name):
        '''
        @brief
        @param basis_name Name of the basis.
        @date 2014
        @author Karl R. Leikanger.
        '''
        self.basis_name = basis_name

    def read_dalton_basisfile(self, filename):
        '''
        @brief Read a basis input file of DALTON format.
        @param filename Name of the inputfile.
        @date 2014
        @author Karl R. Leikanger.
        '''
        # get_words is a "functor". Returning 'words' of the lines 
        # of filename, one line at a time.
        cmarkers = ['!', '$']
        get_words = ReadFile(cmarkers, filename).get_words_of_line
        
        l = -1
        elem = ''

        #iterating over lines in filename 
        for words in iter(get_words, ''):
            if words[0][0] == 'a' or words[0][0] == 'A': 
                if l > -1:
                    self.chem_elems.append(elem)
                l = 0
                elem = self.ChemElemBasis()
                elem.atomic_nr = words[1]
                words = get_words()
            else:
                l += 1
            
            num_exp = int(words[0])
            num_coef = int(words[1])
            # ? = words[2]

            #Store exponents and contraction coeffs to e and c
            e = []; c = []; i = 0
            for i in range(num_exp):
                words = get_words()
                e.append(words[0])
                c.append(words[1:])
                if len(c[-1]) < num_coef:
                    while True:
                        words = get_words()
                        c[-1] = c[-1] + words[:]
                        if len(c[-1]) == num_coef: break

            #set up self.chem_elems.ao
            for i in range(num_coef):
                ao = self.AtomicOrbital()
                indx = []
                for j in range(num_exp):
                    if ( c[j][i] != 0 ): indx.append(j)
                ao.c = ([c[j][i] for j in indx])
                ao.e = ([e[j] for j in indx])
                ao.l = l
                elem.ao.append(ao)

        self.chem_elems.append(elem)

        #sort and check for duplicates
        self.sort_chem_elems()


    def sort_chem_elems(self):   
        '''
        @brief Sort chem elem after atomic numbers and check for duplicates.
        @date 2014
        @author Karl R. Leikanger.
        '''
        
        #set up dict { chem_elems.atomic_nr[i] : i }
        elemset = []
        for elem in self.chem_elems[:]:
            elemset.append(elem.atomic_nr)
        self.chem_elems_indx = dict(zip(elemset, range(len(elemset))))

        #Make sure that the atomic_nr are unique for all chem_elems 
        if len(elemset) > len(set(elemset)):
            print('Error: Check input basisfile format.')
            print('Some chem_elems represented more that once.')
            sys.exit(-1)
        
        #Make sure that chem_elems ars sorted after atomic number 
        if (sorted(elemset) != elemset):
            self.chem_elems = sorted(self.chem_elems, \
                    key=lambda ChemElemBasis: ChemElemBasis.atomic_nr)
            #Dict for sorted list of chem_elems
            self.chem_elems_indx = \
                    dict(zip(sorted(elemset), range(len(elemset))))


    def get_index_of_elem(self, atomic_nr):
        '''
        @brief return the index (in elemset) of the basis for atomic_nr.
        @param atomic_nr Atomic number.
        @date 2014
        @author Karl R. Leikanger.

        '''
        return self.chem_elems_indx.get(atomic_nr)


    def print_cp2k_basisfile(self, filename):
        '''
        @brief return the index (in elemset) of the basis for atomic_nr.
        @param filename Write basis to a file filename.
        @date 2014
        @author Karl R. Leikanger.

        '''

        #open ofilestream

        for elem in self.chem_elems:

            #extract the unique exponents arrays
            exp_unique = []
            for tmp in elem.exp:
                if tmp not in exp_unique:
                    exp_unique.append(e)
        
            #print Element symbol Name of the basis set  Alias names
            #elem_symbol = param.elem_symbols.get(elem.atomic_nr) Basis_name = self.basis_name
            #print nset (repeat the following block of lines nset times)
            #nset=len(exp_unique)
            
            #sort the ao's after the exponents arrays
            aos = sorted(elem.ao, key=lambda AtomicOrbital: AtomicOrbital.exp)

            indx_exp_unique = 0
            indx_aos = 0
            for i in range(len(exp_unique)):

                #collect the aos with the same exponents arrays
                aos_set = []
                while aos[indx_aos].exp == exp_unique[indx_exp_unique]:
                    aos_set.append(aos[indx_aos])
                    indx_aos += 1
                    if (indx_aos==len(aos)): break

                #sort aos after l
                aos_set = sorted(aos_set, key=lambda AtomicOrbital: AtomicOrbital.l)

                indx_exp_unique += 1

                #print set to file
                #print n lmin lmax nexp nshell(lmin) nshell(lmin+1) ... nshell(lmax-1) nshell(lmax)
                #n=1 lmin=aos_set[0].l lmax=aos_set[-1].l [ len(aos_set[k].c) for k in range(len(aos_set)) ]
                #print exponent coeffs(lmin), ..., coeff(lmax)
                #exp_unique[i] aos_set[1].c[i] aos_set[2].c[i] ... aos_set[len(aos_set)-1].c[i] 
                # for i = 1,
                #     i = 2,
                #     ...


        #close ofilestream


basis = Basis('6-311G**')
basis.read_dalton_basisfile('6-311G**')

