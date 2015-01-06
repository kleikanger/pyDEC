'''
@file Classes to handle IO.
'''

import numpy as np
import sys


class FortranIO():
    '''
    @brief Fortran IO superclass.
    @var _int_type Type of intefer numpy.int<xx>
    @var _double_type Type of double numpy.double<xx>
    @var _buf_type numpy.double<xx>
    @var _f filestream
    @date 2014
    @author Karl R. Leikanger.

    Fortran IO is record based, not stream based:
    For unformated IO, Fortran compilers typically write the length of
    the record at the beginning and end of the record. Most but not all
    compilers use four bytes. buf_type must have the same type as the
    buffer.
    '''
    _buf_type = ''
    _types = ''
    _f = ''

    def __init__(self, buf_type, filename, ds):
        '''
        @param buf_type numpy.int<xx>
        @param ds String 'rb', 'wb' for call to open.
        @date 2014
        @author Karl R. Leikanger.
        '''
        self._types = dict()
        self._buf_type = buf_type

        try:
            self._f = open(filename, ds)
        except IOError as e:
            print("Error while trying to open %s", filename)
            print("IO error({0}): {1}".format(e.errno, e.strerror))
            raise SystemExit
        except:
            print('unexcepted error', sys.exc_info()[0])
            raise

    def __del__(self):
        '''
        @date 2014
        @author Karl R. Leikanger.
        '''
        self._f.close()

    def _get_dtype(self, datatype):
        '''
        @param datatype Datatype 'INT', 'DOUBLE', ..., etc.
        @return np datatype
        @date 2014
        @author Karl R. Leikanger.
        '''
        try:
            return self._types[datatype]
        except:
            print('Error: Datatype %s not available?' % datatype)
            raise

    def add_type(self, name, dtype):
        '''
        @brief Add {name: type} to dict _types.
        @param name Name of type.
        @return dtype Data type.
        @date 2014
        @author Karl R. Leikanger.
        '''
        self._types.update({name: dtype})



class ReadFortranBinaryFile(FortranIO):
    '''
    @brief Read fortran (recordbased) binary files.
    @date 2014
    @author Karl R. Leikanger.
    '''

    def __init__(self, buf_type, filename):
        '''
        @date 2014
        @author Karl R. Leikanger.
        '''
        super().__init__(buf_type, filename, 'rb')
        # FortranIO.__init__(self, int_type, double_type, filename, 'rb')
        print('Reading binary file %s.' % filename)

    def read_record(self, datatype, nread):
        '''
        @brief Read a single record from a fortran output file.
        @param datatype Input datatype, 'INT', 'DOUBLE', ..., etc.
        @param nread Number of datatype elements to read.
        @return retvec Numpy array with read data.
        @date 2014
        @author Karl R. Leikanger.
        '''
        buf_type = self._buf_type
        dtype = self._get_dtype(datatype)

        c1 = np.fromfile(self._f, dtype=buf_type, count=1)[0]
        retvec = np.fromfile(self._f, dtype=dtype, count=nread)
        c2 = np.fromfile(self._f, dtype=buf_type, count=1)[0]

        # for record based FORTRAN output, the first and last buffer should have
        # the same value.
        if (c1 != c2):
            print('\nError in ReadFortranBinaryFile input file.\
                  Check output format, datatypes etc.\n')
            raise SystemExit

        return retvec


class WriteFortranBinaryFile(FortranIO):
    '''
    @brief Write to fortran (recordbased) binary files.
    @date 2014
    @author Karl R. Leikanger.
    '''

    def __init__(self, buf_type, filename):
        '''
        @date 2014
        @author Karl R. Leikanger.
        '''
        super().__init__(buf_type, filename, 'wb')
        # FortranIO.__init__(self, int_type, double_type, filename, 'wb')
        print('Writing binary file %s.' % filename)

    def write_record(self, datatype_str, output, buf=True,\
                     print_to_screen=false):
        '''
        @brief Write a single record to a fortran output file.
        @param output Output to write too file.
        @param datatype_str Output datatype, 'INT', 'DOUBLE', ..., etc.
        @param buf Write buffer?
        @param print_to_screen Print to screen?
        @date 2014
        @author Karl R. Leikanger.
        '''
        datatype = self._get_dtype(datatype_str)

        if type(output) == np.ndarray and output.dtype.type != datatype:
            output = output.astype(dtype=datatype)
        else:
            output = np.asarray(output, dtype=datatype)

        if output.ndim == 0:
            output = output.reshape(1)

        if buf:
            buf_type = self._buf_type
            buf = buf_type(output.nbytes)

            self._f.write(buf)
            #buf.tofile(self._f)
            for o in output:
                self._f.write(o)
            # output.tofile(self._f)
            self._f.write(buf)
            #buf.tofile(self._f)

            if print_to_screen:
                print(buf)
                print(output)
                print(buf)
        else:
            for o in output:
                self._f.write(o)
            # output.tofile(self._f)

            if print_to_screen:
                print(output)


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
        for line in iter(self.__f.readline, ''):
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
