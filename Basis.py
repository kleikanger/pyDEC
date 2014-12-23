
'''
@file Contains a class that reads and writes basisfiles in both CP2K and DALTON format.
@brief Read basis file, store basis info and generate basis input to CP2K or DALTON
@date 2014
@author Karl R. Leikanger.

The orbitals are read from file and stored in the following format:
    -chem_elems is an array of __ChemElemBasis objects. Contains atomic number 
    (chem_elems.atomic_nr) and an list of __AtomicOrbital objects (chem_elems.ao), 
    -each of these have a list of coeffs (oa.c), a list of exponents (ao.e) and 
    a AM quantum number (ao.l) stored.  

ex:
filename = '6-311G**'
basis = Basis('basisname')
basis.__read_dalton_basisfile(filename)

'''
import sys
import parameters as param
from IO import ReadFile

class Basis():
    '''
    @brief
    @var chem_elems list of class __ChemElemBasis objects. Basis for given atomic nr. 
    @var chem_elems_indx dictionary {atomic_nr : index in chem_elems }
    @var basis_name Name of basis.
    @date 2014
    @author Karl R. Leikanger.
    '''
    chem_elems = '' 
    chem_elems_indx = '' 
    basis_name = ''
    infile_name = ''

    def __init__(self, basis_name):
        '''
        @date 2014
        @author Karl R. Leikanger.
        '''
        self.basis_name = basis_name
        self.chem_elems = []
    
    class __ChemElemBasis():
        '''
        @brief Container class. Stores basis info for a given element.
        @var atomic_nr Atomic number.
        @var ao List of __AtomicOrbital objects.
        @date 2014
        @author Karl R. Leikanger.
        '''
        atomic_nr = '' 
        ao = ''
        tags_permutations = ''

        def __init__(self):
            self.ao = []
            self.tags_permutations = []

    
    class __AtomicOrbital():
        '''
        @brief Container class. Stores basis info for a given element.
        @var exp Exponents
        @var coef Contraction coefficients
        @var l Quantum number L.
        @var tag To keep track of in which order that the orbitals where read.
        @date 2014
        @author Karl R. Leikanger.
        '''
        e = []
        c = []
        l = ''
        tag = ''

        def __init__(self):
            self.e = []
            self.c = []

    def read_basisfile(self, infile_name, basis_format):
        '''
        @brief Read a basis input file.
        @param infile_name Name of the inputfile.
        @param basis_format Format of inputfile.
        @date 2014
        @author Karl R. Leikanger.
        '''
        self.infile_name = infile_name

        if basis_format in ['DALTON', 'LSDALTON']:
            self.__read_dalton_basisfile()
        else:
            print('Error: Basis input format <%s> not supported.'%basis_format)
            sys.exit(-1)
    
    def print_basisfile(self, outputfile_name, basis_format):
        '''
        @brief Print basis basisfile.
        @param outputfile_name Write basis to a file outputfile_name.
        @param basis_format Format of outputfile.
        @date 2014
        @author Karl R. Leikanger.

        '''

        if basis_format == 'CP2K':
            self.__print_cp2k_basisfile(outputfile_name)
        else:
            print('Error: Output basis format <%s> not supported.'%basis_format)
            sys.exit(-1)

    
    def __read_dalton_basisfile(self):
        '''
        @brief Read a basis input file.
        @date 2014
        @author Karl R. Leikanger.
        '''

        # get_words is a "functor" which returns 'words' of the lines 
        # of infile_name, one line at a time.
        
        print('Reading basis from:', self.infile_name)
        cmarkers = ['!', '$']
        get_words = ReadFile(self.infile_name, cmarkers).get_words_of_line
        
        l = -1
        elem = ''
        aos = ''

        #iterating over lines in infile_name 
        for words in iter(get_words, ''):
            if words[0][0] == 'a' or words[0][0] == 'A': 
                if l > -1:
                    self.chem_elems.append(elem)
                tag = 1
                l = 0
                elem = self.__ChemElemBasis()
                elem.atomic_nr = int(words[1])
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
                e.append(float(words[0]))
                c.append([float(x) for x in words[1:]])
                if len(c[-1]) < num_coef:
                    while True:
                        words = get_words()
                        c[-1] = c[-1] + [float(x) for x in words[:]]
                        if len(c[-1]) == num_coef: break

            #set up self.chem_elems.ao, remove coeffs equal to 0
            for i in range(num_coef):
                ao = self.__AtomicOrbital()
                ao.c = [c[j][i] for j in range(num_exp) if c[j][i] != 0.0]
                ao.e = [e[j] for j in range(num_exp) if c[j][i] != 0.0]
                ao.l = l
                ao.tag = tag; tag += 1
                elem.ao.append(ao)

        self.chem_elems.append(elem)

        #sort and check for duplicates
        self.__sort_chem_elems()


    def __sort_chem_elems(self):   
        '''
        @brief Sort chem elem after atomic numbers and check for duplicates.
        @date 2014
        @author Karl R. Leikanger.
        '''
        #set up list of all Chemical element in the basis
        elemset = []
        for elem in self.chem_elems[:]:
            elemset += [ elem.atomic_nr ]

        #Make sure that the atomic_nr are unique for all chem_elems 
        if len(elemset) > len(set(elemset)):
            print('Error: Check input basisfile format.')
            print('Some chem_elems represented more that once.')
            sys.exit(-1)
        
        #Make sure that chem_elems ars sorted after atomic number 
        if (sorted(elemset) != elemset):
            self.chem_elems = sorted(self.chem_elems, \
                    key=lambda __ChemElemBasis: __ChemElemBasis.atomic_nr)

        #set up dict { chem_elems.atomic_nr[i] : i }
        self.chem_elems_indx = \
                dict(zip(sorted(elemset), range(len(elemset))))

    def __print_cp2k_basisfile(self, filename):
        '''
        @brief Print basis to CP2K basisfile.
        @param filename Write basis to a file filename.
        @date 2014
        @author Karl R. Leikanger.

        '''
        #open ofilestream
        f = open(filename,'w')
        f.write('#\n# CP2K input file.\n')
        f.write('# This file is automatically generated by PROGNAME.\n')
        f.write('# from the DALTON input file %s.\n'%self.infile_name)
        print('Writing basis in CP2K format to file: %s'%filename)

        for elem in self.chem_elems:
            self.__print_elem_to_cp2k_basisfile(elem, f)
        
        f.close()

    def __print_elem_to_cp2k_basisfile(self, elem, f):
        '''
        @brief Print C to CP2K basisfile.
        @param elem __ChemElemBasis element.
        @param f Filestream.
        @date 2014
        @author Karl R. Leikanger.

        Orbitals are written to file in a specific order. 
        They are sorted after:
            - Their exponents arrays.
            - A.M (within each exponents array).

        '''
        #sort the ao's after the exponents arrays
        aos = sorted(elem.ao, key=lambda __AtomicOrbital: __AtomicOrbital.e)

        #extract the unique exponents arrays
        exp_unique = []
        for ao in aos:
            if ao.e not in exp_unique:
                exp_unique.append(ao.e)
            
        ##if any exponents are equal; combine exponents arrays
        #repeat = True
        #while repeat:
        #    repeat = False
        #    for i in range(len(exp_unique)):
        #        for j in range(i+1, range(len(exp_unique))):
        #            si = set(exp_unique[i])
        #            sj = set(exp_unique[j])
        #            if si.intersection(sj) != {}:
        #                #combine exp_unique[i] exp_unique[j] 
        #                i = j = len(exp_unique)
        #                exp_unique[i] = list(si.union(sj))
        #                exp_unique.pop(j)
        #                repeat = True

        #print Element symbol Name of the basis set  Alias names
        elem_symbol = param.elem_symbol_from_a.get(elem.atomic_nr) 
        f.write('#\n# ----------------- \n#\n')
        f.write('%s %s \n'%(elem_symbol, self.basis_name))
        #print nset (repeat the following block of lines nset times)
        f.write('%i\n'%(len(exp_unique)))

        indx_aos = 0
        for i in range(len(exp_unique)):

            #collect the aos with the same exponents arrays
            aos_set = []
            while aos[indx_aos].e == exp_unique[i]:
            #while set(aos[indx_aos].e).intersection(set(exp_unique[i])):
                aos_set.append(aos[indx_aos])
                indx_aos += 1
                if (indx_aos==len(aos)): break #error

            #sort aos after l
            aos_set = sorted(aos_set, key=lambda __AtomicOrbital: __AtomicOrbital.l)

            #collect min and max A.M.
            lmin =  aos_set[0].l
            lmax =  aos_set[-1].l

            #sort coefs after A.M. (in c_l)
            c_l = []
            tmp = []
            ltmp = lmin
            for j in range(lmax-lmin+1):
                if aos_set[j].l != ltmp:
                    c_l.append(tmp)
                    tmp = []
                    ltmp += 1
                tmp.append(aos_set[j].c)
            c_l.append(tmp)

            #get number of elements with each A.M.
            lens = [ len(c_l[k]) for k in range(lmax-lmin+1) ]

            #keep track of the order the orbitals are printed
            elem.tags_permutations += [ x.tag for x in aos_set ]

            #print n lmin lmax nexp nshell(lmin) nshell(lmin+1) ... nshell(lmax)
            f.write('%i %i %i %i '\
                    %(1, aos_set[0].l, aos_set[-1].l, len(exp_unique[i])))
            for ln in lens: 
                f.write('%i '%ln)
            f.write('\n')
            #print exponent coeffs(lmin), ..., coeff(lmax)
            for j in range(len(exp_unique[i])):
                f.write('%f '%(exp_unique[i][j]))
                for l in range(lmax-lmin+1):
                    for k in range(lens[l]):
                        f.write('%f '%c_l[l][k][j])
                f.write('\n')
    
    def __get_index_of_elem(self, atomic_nr):
        '''
        @brief return the index (in elemset) of the basis for atomic_nr.
        @param atomic_nr Atomic number.
        @date 2014
        @author Karl R. Leikanger.

        '''
        return self.chem_elems_indx.get(atomic_nr)


#basis = Basis('STO-6G', 'STO-6G')
#basis.__read_dalton_basisfile()
#basis.__print_cp2k_basisfile('cp2k_STO-6G')

#print(len(basis.chem_elems))
#for e in basis.chem_elems[int(sys.argv[1])].ao:
#    print(e.l)
#    print(e.c)
#    print(e.e)
