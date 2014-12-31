'''
@file Classes to handle IO.
'''

import numpy as np
import sys


class ReadFile():
    '''
    @brief Read lines from file and return as list of strings. Ignore lines that
    @brief starts with cmarkers.
    @var __f Filestream
    @var __cmarkers List of chars. Lines that starts with cmarker are ignored.
    @var __line_nr Wee keep track of the line number.
    @date 2014
    @author Karl R. Leikanger.
    '''
    __f = ''
    __cmarkers = ''
    __line_nr = ''
    __skip_empty = ''

    def __init__(self, filename, cmarkers='', skip_empty=True):
        '''
        @brief Constructor.
        @param __cmarkers As explained above.
        @param __filename Input file.
        @param __skip_empty Skip empty lines?
        @date 2014
        @author Karl R. Leikanger.
        '''
        if type(cmarkers) == str:
            cmarkers = [cmarkers]
        self.__cmarkers = cmarkers
        self.__line_nr = 0
        self.__skip_empty = skip_empty

        try:
            self.__f = open(filename, "r")
        except IOError as e:
            print('Error while trying to open %s' % filename)
            print('IO error({0}): {1}'.format(e.errno, e.strerror))
            raise SystemExit
        except:
            print('unexcepted error', sys.exc_info()[0])
            raise

    def __del__(self):
        '''
        @brief Destructor.
        @date 2014
        @author Karl R. Leikanger.
        '''
        self.__f.close()

    def get_words_of_line(self):
        '''
        @return The next uncommented line as a list of words.
        @date 2014
        @author Karl R. Leikanger.
        '''
        for line in iter(self.f.readline, ''):
            self.__line_nr += 1
            words = line.split()
            if words == []:
                if self.__skip_empty:
                    continue
                else:
                    return ''
            elif words[0][0] in self.__cmarkers:
                continue
            else:
                return words
        return ''

    def get_line_nr(self):
        '''
        @return current line nr of infile.
        @date 2014
        @author Karl R. Leikanger.
        '''
        return self.__line_nr

    def skip_lines(self, n):
        '''
        @brief Skip n lines in file.
        @date 2014
        @param n Nr. of lines to skip.
        @author Karl R. Leikanger.
        '''
        for i in range(n):
            self.__line_nr += 1
            self.__f.readline()


class FortranIO():
    '''
    @var __int_type Type of intefer numpy.int<xx>
    @var __double_type Type of double numpy.double<xx>
    @var __buf_type numpy.double<xx>
    @date 2014
    @author Karl R. Leikanger.

    Fortran IO is record based, not stream based:
    For unformated IO, Fortran compilers typically write the length of
    the record at the beginning and end of the record. Most but not all
    compilers use four bytes. buf_type must have the same type as the
    buffer.
    '''
    __int_type = ''
    __double_type = ''
    __buf_type = ''
    __f = ''

    def __init__(self, int_type, double_type, filename, ds):
        '''
        @param int_type numpy.int<xx>
        @param double_type numpy.double<xx>
        @date 2014
        @author Karl R. Leikanger.
        '''
        self.__int_type = int_type
        self.__double_type = double_type
        self.__buf_type = int_type

        try:
            self.__f = open(filename, ds)
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
        self.__f.close()

    def __get_dtype(self, datatype):
        '''
        @param datatype Datatype 'INT', 'DOUBLE', ..., etc.
        @return np datatype
        @date 2014
        @author Karl R. Leikanger.
        '''
        options = {
            'INT': self.__int_type,
            'DOUBLE': self.__double_type
        }
        return options[datatype]


class ReadFortranBinaryFile(FortranIO):
    '''
    @brief Read fortran (recordbased) binary files.
    @date 2014
    @author Karl R. Leikanger.
    '''

    def __init__(self, int_type, double_type, filename):
        '''
        @date 2014
        @author Karl R. Leikanger.
        '''
        super().__init__(int_type, double_type, filename, 'rb')
        print('Reading binary file %s assuming %iB integer and %iB float.'
              % (filename, int_type(0).nbytes, double_type(0).nbytes))

    def read_record(self, datatype, nread):
        '''
        @brief Read a single record from a fortran output file.
        @param datatype Input datatype, 'INT', 'DOUBLE', ..., etc.
        @param nread Number of datatype elements to read.
        @return retvec Numpy array with read data.
        @date 2014
        @author Karl R. Leikanger.
        '''
        dtype = self.__get_dtype(datatype)
        buf_type = self.__buf_type

        c1 = np.fromfile(self.__f, dtype=buf_type, count=1)[0]
        retvec = np.fromfile(self.__f, dtype=dtype, count=nread)
        c2 = np.fromfile(self.__f, dtype=buf_type, count=1)[0]

        # for record based FORTRAN output, the first and last buffer should have
        # the same value.
        if (c1 != c2):
            # TODO Detailed error message
            print('\nError input file. Check output format, datatypes etc.\n')
            exit(-1)

        return retvec


class WriteFortranBinaryFile(FortranIO):
    '''
    @brief Write to fortran (recordbased) binary files.
    @date 2014
    @author Karl R. Leikanger.
    '''

    def __init__(self, int_type, double_type, filename):
        '''
        @date 2014
        @author Karl R. Leikanger.
        '''
        super().__init__(int_type, double_type, filename, 'wb')
        print('Writing binary file %s assuming %iB integer and %iB float.'
              % (filename, int_type(0).nbytes, double_type(0).nbytes))

    def write_record(self, datatype, output):
        '''
        @brief Write a single record to a fortran output file.
        @param output Output to write too file.
        @param datatype Output datatype, 'INT', 'DOUBLE', ..., etc.
        @date 2014
        @author Karl R. Leikanger.
        '''
        dtype = self.__get_dtype(datatype)

        buf_type = self.__buf_type
        op = dtype(output)
        buf = buf_type(op.nbytes)

        self.__f.write(buf)
        self.__f.write(op)
        self.__f.write(buf)
