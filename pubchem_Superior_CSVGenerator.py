import os


pyScript_directory = os.path.dirname(os.path.abspath("__file__"))

def filePathGenerator(fileDirectoryAndName):
        return os.path.join(pyScript_directory , fileDirectoryAndName)

csvPositive_filePath = filePathGenerator("Debug_CSVfiles/passedSynLigEntries.txt")
csvNegative_filePath = filePathGenerator("Debug_CSVfiles/failedSynLigEntries.txt")
csvPositiveOutput = open(filePathGenerator(csvPositive_filePath),'w')
csvNegativeOutput = open(filePathGenerator(csvNegative_filePath),'w')

inputFileNumber = 25000
multipleFileIterator=0
fileIterator = 25000
linenum = 1

def fileNameToNineDigits(fileNameNumber):
        fileNameLengthBoolean = True
        while fileNameLengthBoolean:
                if len(fileNameNumber) < 9:
                        fileNameNumber = "0" + fileNameNumber
                else:
                        fileNameLengthBoolean = False
        return fileNameNumber

while True:
        print("beginning accepted synonym CSV generation")
        fileIDnum = inputFileNumber + (multipleFileIterator*fileIterator)
        
        inputFileName = "Lig_txt_ready_entries/pubchem_superior-lig-ready_" + fileNameToNineDigits(str(fileIDnum)) + ".txt"
        inputFile = open(filePathGenerator(inputFileName),'r')
        line = inputFile.readline()
        while line != "":
                # ss[0]: "name|syn1|syn2|..."
                # ss[1]: first descriptor; ss[2]: second descriptor;
                # ss[3]: "<SMILES>,x"
                ss = line.split('\t')
                ss[0] = ss[0].strip()
                splitSyns = ss[0].split("|") # splits name from each syn, stored in splitSyns list

                # remove preferred name in first entry- comment out if you want to keep
                splitSyns.pop(0)
                for eachSyn in splitSyns:
                        csvPositiveOutput.write(str(eachSyn)+"\n")
                        csvPositiveOutput.flush()
                #print(str(linenum))
                linenum = linenum + 1
                if linenum%100 == 0:
                        print(str(linenum))
                line=inputFile.readline()
        inputFile.close()
        print("beginning rejected synonym and reason CSV generation")
        inputRejectionFileName = "Debug/pubchem_synonymRejections" + fileNameToNineDigits(str(fileIDnum)) + ".txt"
        inputRejectFile = open(filePathGenerator(inputRejectionFileName),'r')
        lineRejects = inputRejectFile.readline()
        while lineRejects != "":
                # ss[0]: "syn1|syn2|..."
                # ss[1]: "reasonForRejection1|reasonforSyn2|reasonforsyn3|..."
                lineRejects = lineRejects.strip()
                ss = lineRejects.split('\t')
                if len(ss)>1:
                        splitSyns = ss[0].split("|")
                        splitReasons = ss[1].split("|")
                        iterSplitReasons = 0
                        if splitReasons[0].startswith("11: "):
                                combinedOutput = str(eachSyn) + "\t" + str(linenum) + "\t" + "11 with " +str(len(splitSyns))
                        else:
                                try:
                                        for eachSyn in splitSyns:
                                                combinedOutput = str(eachSyn) + "\t" + str(linenum) + "\t" + splitReasons[iterSplitReasons]
                                                iterSplitReasons = iterSplitReasons + 1
                                except:
                                        print("error mismatch-----"+str(len(splitSyns))+"---"+str(len(splitReasons))+"---------------------------------------------")
                        csvNegativeOutput.write(combinedOutput + "\n")
                        csvNegativeOutput.flush()
                        linenum = linenum + 1
                if linenum%100 == 0:
                        print(str(linenum))
                lineRejects = inputRejectFile.readline()  
        multipleFileIterator = multipleFileIterator + 1
        
        print("opening new file number: \t" + str(inputFileNumber + (multipleFileIterator*fileIterator)))
out.close()
