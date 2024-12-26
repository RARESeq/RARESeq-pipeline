import gzip

class FastQRead():
   """A record corresponding to a single read and its associated metadata in a FastQ file"""
   def __init__(self,seqHeader,seqStr,qualHeader,qualStr):
       self.header = seqHeader
       self.seq = seqStr
       self.qualHeader = qualHeader
       self.quals = qualStr
   def toString(self):
       return "%s\n%s\n%s\n%s" %(self.header,self.seq,self.qualHeader,self.quals)

class FastQParser(object):
    """Returns a read-by-read fastQ parser analogous to file.readline()"""
    def __init__(self,filePath,headerSymbols=['@','+']):
        """Returns a read-by-read fastQ parser analogous to file.readline().
        Exmpl: parser.next()
        -OR-
        Its an iterator so you can do:
        for read in parser:
            ... do something with rec ...

        read is a FastQRead object with the fields:
        header,seq,qualHeader,quals
        """
        if filePath.endswith('.gz'):
            self._file = gzip.open(filePath)
        else:
            self._file = open(filePath, 'rU')
        self._currentLineNumber = 0
        self._hdSyms = headerSymbols
        
    def __iter__(self):
        return self
    
    # added by HS 2017-07-06
    def close(self):
      self._file.close()

    def __next__(self):
        return self.next()

    def next(self):
        """Reads in next element, parses, and does minimal verification.
        Returns: FasQRead object with fields: (header,seq,qualHeader,quals)"""
        # ++++ Get Next Four Lines ++++
        elemList = []
        for i in range(4):
            line = self._file.readline()
            #self._currentLineNumber += 1 ## increment file position
            if line:
                elemList.append(line.strip('\n'))
            else:                
                #return
                raise StopIteration 
                #elemList.append(None)
        
    # HS 2017-07-09: commented out below for performance

        ## ++++ Check Lines For Expected Form ++++
        #trues = [bool(x) for x in elemList].count(True)
        #nones = elemList.count(None)
        ## -- Check for acceptable end of file --
        #if nones == 4:
        #   raise StopIteration
        ## -- Make sure we got 4 full lines of data --
        #assert trues == 4,\
        #      "** ERROR: It looks like I encountered a premature EOF or empty line.\n\
        #      Please check FastQ file near line number %s (plus or minus ~4 lines) and try again**" % (self._currentLineNumber)
        ## -- Make sure we are in the correct "register" --
        #assert elemList[0].startswith(self._hdSyms[0]),\
        #      "** ERROR: The 1st line in fastq element does not start with '%s'.\n\
        #      Please check FastQ file near line number %s (plus or minus ~4 lines) and try again**" % (self._hdSyms[0],self._currentLineNumber) 
        #assert elemList[2].startswith(self._hdSyms[1]),\
        #      "** ERROR: The 3rd line in fastq element does not start with '%s'.\n\
        #      Please check FastQ file near line number %s (plus or minus ~4 lines) and try again**" % (self._hdSyms[1],self._currentLineNumber) 
        ## -- Make sure the seq line and qual line have equal lengths --
        #assert len(elemList[1]) == len(elemList[3]), "** ERROR: The length of Sequence data and Quality data of the last record aren't equal.\n\
        #      Please check FastQ file near line number %s (plus or minus ~4 lines) and try again**" % (self._currentLineNumber) 
        
        # ++++ Return fatsQ data as tuple ++++
        return FastQRead(elemList[0],elemList[1],elemList[2],elemList[3])
