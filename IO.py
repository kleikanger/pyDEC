'''
@file Classes to handle IO.
'''

import numpy as np
import sys


class ReadFile():
    '''
    @brief Read lines from file and return as list of strings. Ignore lines that
    @brief starts with cmarkers.
    @var f Filestream
    @var cmarkers List of chars. Lines that starts with cmarker are ignored.
    @var line_nr Wee keep track of the line number.
    @date 2014
    @author Karl R. Leikanger.
    '''
    f = ''
    cmarkers = ''
    line_nr = ''
    skip_empty = ''

    def __init__(self, filename, cmarkers='', skip_empty=True):
        '''
        @brief Constructor.
        @param cmarkers As explained above.
        @param filename Input file.
        @param skip_empty Skip empty lines?
        @date 2014
        @author Karl R. Leikanger.
        '''
        if type(cmarkers) == str:
            cmarkers = [cmarkers]
        self.cmarkers = cmarkers
        self.line_nr = 0
        self.skip_empty = skip_empty

        try:
            self.f = open(filename, "r")
        except IOError as e:
            print('Error while trying to open %s' % filename)
            print('IO error({0}): {1}'.format(e.errno, e.strerror))
            sys.exit(-1)
        except:
            print('unexcepted error', sys.exc_info()[0])
            raise

    def __del__(self):
        '''
        @brief Destructor.
        @date 2014
        @author Karl R. Leikanger.
        '''
        self.f.close()

    def get_words_of_line(self):
        '''
        @brief Returns the next uncommented line as a list of words.
        @date 2014
        @author Karl R. Leikanger.
        '''
        for line in iter(self.f.readline, ''):
            self.line_nr += 1
            words = line.split()
            if words == []:
                if self.skip_empty:
                    continue
                else:
                    return ''
            elif words[0][0] in self.cmarkers:
                continue
            else:
                return words
        return ''

    def get_line_nr(self):
        return self.line_nr

    def skip_lines(self, n):
        '''
        @brief Skip n lines in file.
        @date 2014
        @param n Nr. of lines to skip.
        @author Karl R. Leikanger.
        '''
        for i in range(n):
            self.line_nr += 1
            self.f.readline()


class ReadFortranBinaryFile():
    '''
    @brief Read fortran (recordbased) binary files.
    @var int_type Type of intefer numpy.int<xx>
    @var double_type Type of double numpy.double<xx>
    @var buf_type numpy.double<xx>
    @date 2014
    @author Karl R. Leikanger.

    Fortran IO is record based, not stream based:
    For unformated IO, Fortran compilers typically write the length of
    the record at the beginning and end of the record. Most but not all
    compilers use four bytes. buf_type must have the same type as the
    buffer.
    '''

    int_type = ''
    double_type = ''
    buf_type = ''
    f = ''

    def __init__(self, int_type, double_type, filename):
        '''
        @param int_type numpy.int<xx>
        @param double_type numpy.double<xx>
        @date 2014
        @author Karl R. Leikanger.
        '''
        self.int_type = int_type
        self.double_type = double_type
        self.buf_type = int_type

        print('Reading binary file %s assuming %iB integer and %iB float.'
              % (filename, int_type(0).nbytes, double_type(0).nbytes))

        try:
            self.f = open(filename, "rb")
        except IOError as e:
            print("Error while trying to open %s", filename)
            print("IO error({0}): {1}".format(e.errno, e.strerror))
            sys.exit(-1)
        except:
            print('unexcepted error', sys.exc_info()[0])
            raise

    def __del__(self):
        '''
        @date 2014
        @author Karl R. Leikanger.
        '''
        self.f.close()

    def read_record(self, datatype, nread):
        '''
        @brief Read a single record from a fortran output file.
        @param f Inputstream.
        @param datatype Numpy datatype.
        @param nread Number of datatype elements to read.
        @return retvec Numpy array with read data.
        @date 2014
        @author Karl R. Leikanger.

        '''
        if datatype == 'INT':
            dtype = self.int_type
        elif datatype == 'DOUBLE':
            dtype = self.double_type

        buf_type = self.buf_type
        f = self.f

        c1 = np.fromfile(f, dtype=buf_type, count=1)[0]
        retvec = np.fromfile(self.f, dtype=dtype, count=nread)
        c2 = np.fromfile(self.f, dtype=buf_type, count=1)[0]

        # for record based FORTRAN output, the first and last buffer should have
        # the same value.
        if (c1 != c2):
            # TODO Detailed error message
            print('\nError in input file. Check output format, datatypes etc.\n')
            exit(-1)

        return retvec
