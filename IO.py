
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

    def __init__(self, cmarkers, filename):
        '''
        @brief Constructor.
        @param cmarkers As explained above.
        @param filename Input file.
        @date 2014
        @author Karl R. Leikanger.
        '''
        self.cmarkers = cmarkers
        self.line_nr = 0
    
        print('Reading basis from:', filename)
        try:
            self.f = open(filename, "r")
        except IOError as e:
            print("Error while trying to open %s", self.filename)
            print("IO error({0}): {1}".format(e.errno, e.strerror))
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
            if words == []: continue
            if words[0][0] in self.cmarkers: continue
            else: return words
        return ''

