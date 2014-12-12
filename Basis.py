
'''
@brief Read basis file, store basis info and generate basis input to CP2K or DALTON
@date 2014
@author Karl R. Leikanger.

Basis set format is equivalent to the one in CP2K.
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

class Basis():
    '''
    @brief
    @var elements list of class ElemBasis objects. Basis for given atomic nr. 
    @var elems_dict dictionary {atomic_nr : index in elements }
    @date 2014
    @author Karl R. Leikanger.
    '''

    elements = [] 
    elems_dict = '' 
    
    class ElemBasis()
        '''
        @brief Container class. Stores basis info for a given element.
        @var a Atomic number
        @var nset Number of exponent sets
        @var n Principle quantum number
        @var lmax Maximum angular momentum quantum number l
        @var lmin Minimum angular momentum quantum number l
        @var nshell Number of shells for angular momentum quantum number l
        @var e Exponent
        @var c Contraction coefficients
        @date 2014
        @author Karl R. Leikanger.
        '''
        a = '' 
        num_aos = ''
        num_coef = ''
        e = []
        c = []
        s = []
        l = []
        n = '' 


        #nset = ''
        #lmax = ''
        #lmin = ''
        #nshell = []
    
    def __init__(self, fname):
        '''
        @brief
        @param
        @date 2014
        @author Karl R. Leikanger.
        '''

    def read_dalton_basisfile():
        '''
        @brief
        @param
        @date 2014
        @author Karl R. Leikanger.
        '''

        try:
            f = open(self.filename, "r")
        except IOError as e:
            print("Error while trying to open %s", self.filename)
            print("IO error({0}): {1}".format(e.errno, e.strerror))
        except:
            print('unexcepted error', sys.exc_info()[0])
            raise

        l = 0
        for line in f.readline():
            if line[0] == '$' continue
            elif line[0] == 'a' or line[0] == 'A':
                
                words = line.split()
                elem = ElemBasis()
                elem.a = word[1]
                elem.l = l
                l += 1

                words = f.readline().split()
                elem.num_aos = words[0]
                elem.num_coef = words[1]
                # ? = words[2]

                for i in range(elem.num_aos):
                    words = f.readline().split()
                    elem.e.append(words[0])
                    elem.c.append(words[1:elem.num_coef])

                    #TODO
                    #elem.n=append

                self.elements.append(elem)
            else:
                print('Error: Check input basisfile format.')
                sys.exit(-1)
        
        elemset = []
        for elem in self.elem[:]:
            elemset.append(elem.a)
        if len(elemset) > len(set(elemset)):
            print('Error: Check input basisfile format.')
            print('Some elements represented more that once.')
            sys.exit(-1)
 
        self.elems_dict = dict( zip(elemset, range( len(elemset) )) )

    def get_index_of_elem(self, atomic_nr):
        '''
        @brief return the index (in elemset) of the basis for atomic_nr.
        @param atomic_nr Atomic number.
        @date 2014
        @author Karl R. Leikanger.

        '''
        return self.elems_dict.get(atomic_nr)
        

            






        


