'''
@file Contains a class that reads/writes basisfiles in both CP2K/DALTON format.
@brief Read basis file, store basis info and generate basis input to CP2K or
DALTON
@date 2014
@author Karl R. Leikanger.

The orbitals are read from file and stored in the following format:
    -chem_elems is an array of __ChemElemBasis objects. Contains atomic number
    (chem_elems.atomic_nr) and an list of __AtomicOrbital objects
    (chem_elems.ao),
    -each of these have a list of coeffs (oa.c), a list of exponents (ao.e) and
    a AM quantum number (ao.l) stored.
'''
import Parameters as param
from IO import ReadFile
import numpy as np


class BasisSet():
    '''
    @brief
    @var chem_elems list of class __ChemElemBasis objects. One for each element.
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

        def __init__(self):
            self.ao = []

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

        options = {
            'dalton': self.__read_dalton_basisfile,
            'lsdalton': self.__read_dalton_basisfile
        }
        try:
            options[basis_format.lower()]()
        except:
            print('Error: Basis input format <%s> not supported?'
                  % basis_format)
            raise  # SystemExit

    def print_basisfile(self, basis_format, f):
        '''
        @brief Print basis basisfile.
        @param f Output filestream.
        @param basis_format Format of outputfile.
        @date 2014
        @author Karl R. Leikanger.
        '''
        options = {
            'cp2k': self.__print_cp2k_basisfile,
            'dalton': self.__print_dalton_basisfile,
            'lsdalton': self.__print_dalton_basisfile
        }
        try:
            options[basis_format.lower()](f)
        except:
            print('Error: Output basis format <%s> not supported?'
                  % basis_format)
            raise  # SystemExit

    # TODO not needed? delete?
    def get_mo_transformation(self, codeformat_from, codeformat_to, atoms):
        '''
        @brief Set up the transformation of MO's from codeformat_to to\
            codeformat_from.
        The permutation of the mo's is then simple:
            - mo = [mo[i] for i in permutations]
        @return The output is a list with permutations.
        @param codeformat_from String object: ('LSDALTON', 'CP2K', ...)
        @param codeformat_to String object: ('LSDALTON', 'CP2K', ...)
        @param atoms Elements from inputfile type System.__Atom
        @date 2014
        @author Karl R. Leikanger.
        '''
        options = {
            'dalton': self.__get_aoorder_dalton,
            'lsdalton': self.__get_aoorder_dalton,
            'cp2k': self.__get_aoorder_cp2k
        }
        try:
            perm_from = options[codeformat_from.lower()](atoms)
            perm_to = options[codeformat_to.lower()](atoms)
        except:
            print('Error in set_up_ao_order: codeformat <%s> or <%s> not \
                  supported?.' % (codeformat_to, codeformat_from))
            raise

        # duplicate elements or different length of perm_to and perm_from?
        lens = [len(set(perm_to)), len(perm_to),
                len(set(perm_from)), len(perm_from)]
        if min(lens) != max(lens):
            print('Error in set_up_ao_order: Something wrong with the \
                  permutation arrays')
            raise SystemExit

        permutations = []
        for x in perm_to:
            permutations.append(perm_from.index(x))
        return permutations

    def get_ao_order(self, codeformat, atoms):
        '''
        @brief Get order of AO's relative to CODENAME storage of atoms.
        @return The output is a permutation list.
        @param codeformat String object: ('LSDALTON', 'CP2K', ...)
        @param atoms Elements from inputfile type System.__Atom
        @date 2014
        @author Karl R. Leikanger.
        '''
        options = {
            'dalton': self.__get_aoorder_dalton,
            'lsdalton': self.__get_aoorder_dalton,
            'cp2k': self.__get_aoorder_cp2k
        }
        try:
            return options[codeformat.lower()](atoms)
        except:
            print('Error BasisSet:get_ao_order: codeformat <%s> not \
                  supported?.' % (codeformat))
            raise

    def get_index_of_elem(self, atomic_nr):
        '''
        @brief return the index (in elemset) of the basis for atomic_nr.
        @param atomic_nr Atomic number.
        @date 2014
        @author Karl R. Leikanger.

        '''
        return self.chem_elems_indx.get(atomic_nr)

    def get_chem_elem(self, atomic_nr):
        '''
        @brief Find and return the correct __ChemElemBasis object.
        @param atomic_nr Atomic number.
        @return the __ChemElemBasis object for elem with object.a = atomic_nr.
        @date 2014
        @author Karl R. Leikanger.

        '''
        indx = self.chem_elems_indx.get(atomic_nr)
        return self.chem_elems[indx]

    def __get_ao_from_tag(self, elem, tag):
        '''
        @brief Return __AtomicOrbital object with obj.tag==tag
        @param elem Chemical element.
        @param tag Tag of ao.
        @return __AtomicOrbital ao with ao.tag==tag
        @date 2014
        @author Karl R. Leikanger.
        '''
        # TODO set un dict { tag : ao (ao.tag==tag)}
        for ao in elem.ao:
            if ao.tag == tag:
                return ao
        # TODO setup tagmax with error msg
        print('error: ao with tag %s does not excist' % tag)
        raise SystemExit

    def __sort_chem_elems(self):
        '''
        @brief Sort chem elem after atomic numbers and check for duplicates.
        @date 2014
        @author Karl R. Leikanger.
        '''
        # set up list of all Chemical element in the basis
        elemset = []
        for elem in self.chem_elems[:]:
            elemset += [elem.atomic_nr]

        # Make sure that the atomic_nr are unique for all chem_elems
        if len(elemset) > len(set(elemset)):
            print('Error: Check input basisfile format.')
            print('Some chem_elems represented more that once.')
            raise SystemExit

        # Make sure that chem_elems ars sorted after atomic number
        if (sorted(elemset) != elemset):
            self.chem_elems = sorted(
                self.chem_elems,
                key=lambda __ChemElemBasis: __ChemElemBasis.atomic_nr)

        # set up dict { chem_elems.atomic_nr[i] : i }
        self.chem_elems_indx = \
            dict(zip(sorted(elemset), range(len(elemset))))

    def __read_dalton_basisfile(self):
        '''
        @brief Read a basis input file.
        @date 2014
        @author Karl R. Leikanger.
        '''
        print('Reading basis from:', self.infile_name)
        cmarkers = ['!', '$']
        get_words = ReadFile(self.infile_name, cmarkers).get_words_of_line

        l = -1
        elem = ''

        # iterating over lines in infile_name
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

            # Store exponents and contraction coeffs to e and c
            e = []
            c = []
            i = 0
            for i in range(num_exp):
                words = get_words()
                e.append(float(words[0]))
                c.append([float(x) for x in words[1:]])
                if len(c[-1]) < num_coef:
                    while True:
                        words = get_words()
                        c[-1] = c[-1] + [float(x) for x in words[:]]
                        if len(c[-1]) == num_coef:
                            break

            # set up self.chem_elems.ao, remove coeffs equal to 0
            for i in range(num_coef):
                ao = self.__AtomicOrbital()
                ao.c = [c[j][i] for j in range(num_exp) if c[j][i] != 0.0]
                ao.e = [e[j] for j in range(num_exp) if c[j][i] != 0.0]
                ao.l = l
                ao.tag = tag
                tag += 1
                elem.ao.append(ao)

        self.chem_elems.append(elem)

        # sort and check for duplicates
        self.__sort_chem_elems()

    def __print_dalton_basisfile(self, f):
        '''
        @brief Print basis to DALTON basisfile.
        @param f Output filestream
        @date 2014
        @author Karl R. Leikanger.

        '''
        for elem in self.chem_elems:
            self.__print_elem_to_dalton_basisfile(elem, f)

    def __print_elem_to_dalton_basisfile(self, elem, f):
        '''
        @brief Print __ChemElemBasis object to DALTON basisfile.
        @param elem __ChemElemBasis element.
        @param f Filestream.
        @date 2014
        @author Karl R. Leikanger.
        '''
        # print Element symbol Name of the basis set
        elem_symbol = param.elem_symbol_from_a.get(elem.atomic_nr)
        f.write('$\n$\n$ --- Elem = %s, Basis = %s --- \n'
                % (elem_symbol, self.basis_name))
        f.write('A %i\n' % elem.atomic_nr)

        # get AO's tags sorted in dalton order
        tags = self.__get_order_of_dalton_aos_tags(elem)
        # sort aos acc. to tags
        # indices = [x - min(tags[:]) for x in tags]
        # aos = [elem.ao[i] for i in indices]
        aos = [self.__get_ao_from_tag(elem, tag) for tag in tags]

        lsym = ['S', 'P', 'D', 'F', 'G', 'H', 'I', 'J', 'K', 'L', 'M']
        for l in set([ao.l for ao in aos]):
            f.write('$\n$ %s - type functions:\n' % lsym[l])

            # extract ao's all with the correct AM
            aos_set = [ao for ao in aos if ao.l == l]
            # extract all 'unique' exponents arrays
            exp_unique = []
            for ao in aos_set:
                if ao.e not in exp_unique:
                    exp_unique.append(ao.e)
            # get total nr of exponents
            total_n_exp = 0
            for eu in exp_unique:
                total_n_exp += len(eu)

            # XXX what is the meaning os last integer?? Ask S.R.
            f.write('%i %i %i \n' % (total_n_exp, len(aos_set), 0))

            # organize and print exponents + coefficients
            zlists = np.zeros([total_n_exp, len(aos_set)+1], dtype=float)
            b = 0
            for eu in exp_unique:
                leu = len(eu)
                zlists[b:b+leu, 0] = eu
                b += leu
            k = 0
            b = 0
            for eu in exp_unique:
                leu = len(eu)
                for ao in aos_set:
                    if ao.e == eu:
                        zlists[b:b+leu, k+1] = aos_set[k].c
                        k += 1
                b += leu
            for i in range(zlists.shape[0]):
                for j in range(zlists.shape[1]):
                    ostr = str(zlists[i, j]).ljust(14)
                    f.write('%s' % ostr)
                f.write('\n')

    def __print_cp2k_basisfile(self, f):
        '''
        @brief Print basis to CP2K basisfile.
        @param f Output filestream.
        @date 2014
        @author Karl R. Leikanger.

        '''
        for elem in self.chem_elems:
            self.__print_elem_to_cp2k_basisfile(elem, f)

    def __print_elem_to_cp2k_basisfile(self, elem, f):
        '''
        @brief Print __ChemElemBasis object to CP2K basisfile.
        @param elem __ChemElemBasis element.
        @param f Filestream.
        @date 2014
        @author Karl R. Leikanger.

        Orbitals are written to file in a specific order.
        They are sorted after:
            - Their exponents arrays.
            - A.M (within each exponents array).
        '''
        # tags of the
        tags, n_unique_exp = self.__get_order_of_CP2K_aos_tags(elem)

        # print Element symbol Name of the basis set  Alias names
        elem_symbol = param.elem_symbol_from_a.get(elem.atomic_nr)
        f.write('#\n# ----------------- \n#\n')
        f.write('%s %s \n' % (elem_symbol, self.basis_name))

        # print nset (repeat the following block of lines nset times)
        f.write('%i\n' % n_unique_exp)

        # sort aos after tags
        indices = [x - min(tags[:]) for x in tags]
        aos = [elem.ao[i] for i in indices]

        # iterate over the aos with the same exponents
        i = 0
        for _ in range(n_unique_exp):
            # collect all ao's with equal exponents in aos_set
            e = aos[i].e
            aos_set = []
            for ao in aos[i:]:
                if ao.e != e:
                    break
                aos_set.append(ao)
                i += 1

            # collect min and max A.M.
            lmin = aos_set[0].l
            lmax = aos_set[-1].l

            # sort coefs after A.M. (in c_l)
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

            # get number of elements with each A.M.
            lens = [len(c_l[k]) for k in range(lmax-lmin+1)]

            # print n lmin lmax nexp nshell(lmin) nshell(lmin+1) .. nshell(lmax)
            f.write('%i %i %i %i '
                    % (1, aos_set[0].l, aos_set[-1].l, len(e)))
            for ln in lens:
                f.write('%i ' % ln)
            f.write('\n')

            # print exponent coeffs(lmin), ..., coeff(lmax)
            for j in range(len(e)):
                ostr = str(e[j]).ljust(14)
                f.write('%s ' % ostr)
                for l in range(lmax-lmin+1):
                    for k in range(lens[l]):
                        ostr = str(c_l[l][k][j]).ljust(14)
                        f.write('%s ' % ostr)
                f.write('\n')

    def __get_aoorder_dalton(self, atoms):
        '''
        @brief Set up the order of the a.o.'s as it is stored by the code CP2K.
        @return The output is a list of permutations to pyQCSlink format.
        @param atoms Elements from inputfile type System.__Atom.\
            The atoms in the system in the order they are read by the code.
        @date 2014
        @author Karl R. Leikanger.
        '''
        # How is the orbitals ordered for a given l?
        # l=1:px,py,px l=2:pxy,pxz,pyz,p(xx-yy),pzz etc..
        m_order = {  # l : [orbitals, ...]
            0: [0],
            1: [0, 1, 2],  # x,z,y ??
            # 2 : ? TODO FIXME find the order, see cp2k source.
        }

        # extend tags with px,py,pz,...
        # and append tagslists of different atoms.
        ao_list = []
        listtot = 0
        for atom in atoms:
            # add func for several basis sets.
            elem = self.get_chem_elem(atom.a)  # atom.a = Atomic Number
            tags = self.__get_order_of_dalton_aos_tags(elem)
            order = tags.copy()
            j = -1
            for i in range(len(tags)):
                j += 1
                ao = self.__get_ao_from_tag(elem, tags[i])
                tmp = [(tags[i] + x) for x in m_order[ao.l]]
                for k in range(len(order)):
                    if order[k] > order[j]:
                        order[k] += len(tmp)-1
                order[j] = tmp[0]
                for k in range(len(tmp)-1):
                    order.insert(i+1, tmp[-1-k])
                    j += 1
            ao_list += [(x + listtot) for x in order]
            listtot += max(order)
        return ao_list

    def __get_order_of_dalton_aos_tags(self, elem):
        '''
        @brief Set up the order of the a.o.'s as it is read by the code DALTON.
        @param atoms Elements from inputfile type System.__Atom
        @return Tags of the aos in the stored order of CP2K\
            number of unique exponents arrays.
        @date 2014
        @author Karl R. Leikanger.
        '''
        # sort the ao's after the exponents arrays
        aos = sorted(elem.ao, key=lambda __AtomicOrbital: __AtomicOrbital.e,
                     reverse=True)

        # and then after l
        tags = []
        lmax = max([ao.l for ao in aos])
        for l in range(lmax+1):
            for ao in aos:
                if ao.l == l:
                    tags.append(ao.tag)

        return tags

    def __get_aoorder_cp2k(self, atoms):
        '''
        @brief Set up the order of the a.o.'s as it is stored by the code CP2K.
        @return The output is a list of permutations to pyQCSlink order.
        @param atoms Elements from inputfile type System.__Atom.\
            The atoms in the system in the order they are read by the code.
        @date 2014
        @author Karl R. Leikanger.
        '''
        # see a cp2k/cp2k/src/atoms_imput.F

        # How is the orbitals ordered for a given l
        # l=1:px,py,px l=2:pxy,pxz,pyz,p(xx-yy),pzz etc..
        m_order = {  # l : [orbitals, ...]
            0: [0],
            1: [1, 2, 0],  # x,z,y
            # 2 : ? TODO FIXME find the order, see cp2k source.
        }

        # extend tags with px,py,pz,...
        # and append tagslists of different atoms.
        ao_list = []
        listtot = 0
        for atom in atoms:
            # add func for several basis sets.
            elem = self.get_chem_elem(atom.a)  # atom.a = Atomic Number
            tags, n = self.__get_order_of_CP2K_aos_tags(elem)

            order = tags.copy()
            j = -1

            for i in range(len(tags)):
                j += 1
                ao = self.__get_ao_from_tag(elem, tags[i])
                tmp = [(tags[i] + x) for x in m_order[ao.l]]
                for k in range(len(order)):
                    if order[k] > order[j]:
                        order[k] += len(tmp)-1
                order[j] = tmp[0]
                for k in range(len(tmp)-1):
                    order.insert(i+1, tmp[-1-k])
                    j += 1
            ao_list += [(x + listtot) for x in order]
            listtot += max(order)
        return ao_list

    def __get_order_of_CP2K_aos_tags(self, elem):
        '''
        @brief Set up the order of the a.o.'s as it is read by the code CP2K.
        @param atoms Elements from inputfile type System.__Atom
        @return Tags of the aos in the stored order of CP2K\
            number of unique exponents arrays.
        @date 2014
        @author Karl R. Leikanger.
        '''
        # sort the ao's after the exponents arrays
        aos = sorted(elem.ao, key=lambda __AtomicOrbital: __AtomicOrbital.e)

        # extract the unique exponents arrays
        exp_unique = []
        for ao in aos:
            if ao.e not in exp_unique:
                exp_unique.append(ao.e)

        # XXX if any exponents are equal; combine exponents arrays ?

        indx_aos = 0
        tags = []
        for i in range(len(exp_unique)):

            # collect the aos with the same exponents arrays
            aos_set = []
            while aos[indx_aos].e == exp_unique[i]:
                # while set(aos[indx_aos].e).intersection(set(exp_unique[i])):
                aos_set.append(aos[indx_aos])
                indx_aos += 1
                if (indx_aos == len(aos)):
                    break  # error

            # sort aos after l
            aos_set =\
                sorted(aos_set, key=lambda __AtomicOrbital: __AtomicOrbital.l)

            # the CP2K order of the a.o.'s
            tags += [x.tag for x in aos_set]

        return tags, len(exp_unique)
