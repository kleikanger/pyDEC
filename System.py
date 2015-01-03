'''
@file Store info about the physical system.
@date 2014
@author Karl R. Leikanger
'''

class System():
   '''
   @brief Stores information about reference cell atoms and lattice vectors.
   @var latvec Lattice vectors
   @var atoms Atoms object list.
   @var atombasis Different basis for different elements?
   @var angstrom Input in units of Angstrom instead of Bohr?
   @date 2014
   @author Karl R. Leikanger
   '''
   latvec = ''
   atoms = ''
   atombasis = False
   angstrom = False

   def __init__(self):
      '''
      @param
      @date 2014
      @author Karl R. Leikanger
      '''
      self.latvec = [[0, 0, 0], [0, 0, 0], [0, 0, 0]]
      self.atoms = []

   class __Atom():
      '''
      @brief Atoms storage container
      @var basis Basis set.
      @var a Atomic number.
      @var ghost Is ghost?
      @var position Position of nucleus.
      @date 2014
      @author Karl R. Leikanger
      '''
      atom_tag = ''
      basisname = ''
      a = ''
      position = ''
      ghost = ''
      charge = ''

   def read_system_info(self, filename, fformat):
      '''
      @brief Read dalton molfile and store info.
      @param filename Filename of dalton .mol file.
      @param fformat Input format.
      @date 2014
      @author Karl R. Leikanger
      '''
      if fformat in ['DALTON', 'LSDALTON']:
          self.__read_dalton_molfile(filename)
      else:
        print('Error: System file format <%s> not implemented.' % fformat)

   def __read_dalton_molfile(self, filename):
      '''
      @brief Read dalton molfile and store info.
      @param filename Filename of dalton .mol file.
      @date 2014
      @author Karl R. Leikanger
      '''
      print('Reading DALTON .mol file: %s' % filename)
      fs = ReadFile(filename, skip_empty=False)
      get_words = fs.get_words_of_line
      line_nr = fs.get_line_nr
      skip_lines = fs.skip_lines

      basis = ''

      # first line : BASIS or ATOMBASIS
      word = get_words()[0]
      options = { 'basis' : False, 'atombasis' : True}
      self.atombasis = options[word.lower()]
      if not self.atombasis:
         word = get_words()[0]
         basis = word

      # next two lines are reserved for comments
      skip_lines(2)

      # next line: info about total molecule
      self.charge = 0
      words = get_words()
      for word in words:
         tmp = word.split(sep='=')
         if tmp[0].lower() == 'angstrom':
            self.angstrom = True
         elif tmp[0].lower() == 'charge':
            self.charge = float(tmp[1])
         elif tmp[0].lower() == 'atomtypes':
            self.atomtypes = int(tmp[1])
         else:
            print('Error: Invalid keyword <%s> in line %i of file: %s.'
                  % (tmp[0], line_nr(), filename))
            sys.exit(-1)

      # Next line: info about next atoms or the latticevectors
      tag = 0
      for words in iter(get_words, ''):
         if words == []:
             continue
         natoms = 0
         for word in words:
            tmp = word.split(sep='=')
            if tmp[0].lower() == 'charge':
               charge = float(tmp[1])
            elif tmp[0].lower() == 'atoms':
               natoms = int(tmp[1])
            elif tmp[0].lower() == 'basis':
               if not self.atombasis:
                  print('Error: Atombasis is not set in file: %s => Basis is\
                        not a valid keyword in line %i' % (filename, line_nr()))
                  sys.exit(-1)
               basis = tmp[1]
            elif tmp[0].lower() in ['a1', 'a2', 'a3']:
                indx = int(tmp[0][1]) - 1
                try:
                    self.latvec[indx] = [float(words[i]) for i in [-4, -3, -2]]
                except:
                    print('Error: Input from file: %s, line %s.\n'
                          % (filename, line_nr())
                          + 'Cannot read lattice vectors.')
                    raise
                break
            else:
               print('Error: Invalid keyword <%s> in line %i of file: %s.'
                     % (tmp[0], line_nr(), filename))
               sys.exit(-1)

            # check that all param are set: basis, natom, ...

            # read atoms symbol xpos ypos zpos
            for i in range(natoms):

               atom = self.__Atom()
               atom.basisname = basis
               atom.charge = charge

               tag += 1
               atom.tag = tag

               words = get_words()
               atom.a = parameters.elem_a_from_symbol.get(words[0])
               if atom.a is None:
                  print('Error: Input from file: %s, line %s. Element %s not\
                        in list of symbols.' % (filename, line_nr(), words[0]))
                  exit(0)
               try:
                  atom.position = [float(x) for x in words[1:]]
               except:
                  print('Error: Input from file: %s, line %s.\
                        Cannot read atom coords.' % (filename, line_nr()))
                  raise
               atom.ghost = False

               self.atoms.append(atom)

      for x in self.atoms:
           print(x.basisname)
           print(x.a)
           print(x.position)
           print(x.ghost)
           print(x.charge)
           print('\n')
      for x in self.latvec:
         print(x)

#system = System()
#system.read_dalton_molfile('H2.mol')
