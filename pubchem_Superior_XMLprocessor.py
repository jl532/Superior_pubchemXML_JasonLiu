import xml.etree.ElementTree as etree
from ChemScript17 import *
import re
import time
import os

pyScript_directory = os.path.dirname(os.path.abspath("__file__"))
synonymInput_path = os.path.join(pyScript_directory , "Synonym_File/CID-Synonym-filtered.txt")
filter_begin_path = os.path.join(pyScript_directory , "Filter_Settings_and_Criterion/filter_beginnings.txt")
filter_middle_path = os.path.join(pyScript_directory , "Filter_Settings_and_Criterion/filter_middles.txt")
filter_end_path = os.path.join(pyScript_directory , "Filter_Settings_and_Criterion/filter_endings.txt")

pubChemSyns = open(synonymInput_path, 'r')
filterBegins = open(filter_begin_path,'r')
filterMids = open(filter_middle_path,'r')
filterEnds = open(filter_end_path,'r')

synStatsFile = open("pubchem_Superior_generatedSynonymStatistics.txt",'w')

## filter parameter file processing

filterBeginsReadout = filterBegins.readline()
filterMidsReadout = filterMids.readline()
filterEndsReadout = filterEnds.readline()

def filterFiler(inputReadout):
        filter_strings=[]
        for each in inputReadout.split(","):
                filter_strings.append(each.strip("\'"))
        return filter_strings

filter_beginningStrings = filterFiler(filterBeginsReadout)
filter_middleStrings = filterFiler(filterMidsReadout)
filter_endStrings = filterFiler(filterEndsReadout)

catalogNumRgx1=re.compile("[0-9]{3}") # all capitals/numbers, possible space or hyphen, and all capitals/numbers
catalogNumRgx2=re.compile("[A-Z]{14}-[A-Z]{10}-[A-Z]") #14 char - 10 char - 1 char all caps
casNumRgx=re.compile("[0-9]+-[0-9]+-[0-9]+") # Cas number identifier-- 100% accurate
casNumRgxAlpha=re.compile("[A-Za-z]")

def synFilter(inputSyn):
        outputLogic = [True, "default cases passes", 0]
        if inputSyn == "" or len(inputSyn.strip(" ")) < 3:
                outputLogic = [False , "empty, or too few characters for a useful synonym entry", 1 ]
        else:
                if catalogNumRgx1.search(inputSyn):
                        outputLogic = [False , "fails 3-consecutive-digit test: likely catalog number", 2 ]
                if catalogNumRgx2.search(inputSyn):
                        outputLogic = [False , "fails 14-10-1 digit test: likely catalog string", 3 ]
                if casNumRgx.match(inputSyn):
                        outputLogic = [True , "passes CAS # characteristics", 4 ]
                        if casNumRgxAlpha.search(inputSyn):
                                outputLogic = [False , "fails CAS test: extra characters" , 5 ]
                if inputSyn.startswith(tuple(filter_beginningStrings)):
                        outputLogic = [False , "fails beginning char test: likely catalog number", 6 ]
                for eachRem in filter_middleStrings:
                        if eachRem in inputSyn:
                                outputLogic = [False , "fails middle char test: likely catalog number", 7 ]
                if inputSyn.endswith(tuple(filter_endStrings)):
                        outputLogic = [False , "fails ending char test: likely catalog number", 8 ]
                if len(inputSyn) < 70:
                        if StructureData.LoadData(inputSyn,"name") != None:
                                outputLogic = [False , "fails: ChemDraw recognizes entry, therefore is left out", 9 ]
                else:
                        outputLogic = [False , "fails: input has too many characters" , 10 ]
        return outputLogic

    ## add writing to error debug file to keep track of rejected synonyms
    ## add debug file with statistics about the synonyms accepted/rejected

def fileNameToNineDigits(fileNameNumber):
        fileNameLengthBoolean = True
        while fileNameLengthBoolean:
                if len(fileNameNumber) < 9:
                        fileNameNumber = "0" + fileNameNumber
                else:
                        fileNameLengthBoolean = False
        return fileNameNumber
        

cid = 0
multipleFileIterator=0
fileIterator = 25000
while True:
        
        fileID_Start = str(1 + multipleFileIterator*fileIterator)
        fileID_End = str(25000 + multipleFileIterator*fileIterator)
        fileID_Start_serialized = fileNameToNineDigits(fileID_Start)
        fileID_End_serialized = fileNameToNineDigits(fileID_End)
        
        xmlFileInput = "XML_CID_Data/Compound_" + fileID_Start_serialized + "_" + fileID_End_serialized + ".xml"
        outputFileName = 'Lig_txt_ready_entries/pubchem_superior-lig-ready_' + fileID_End_serialized + '.txt'
        debugFileName = "Debug/pubchem_failureToPassFilter" + fileID_End_serialized + ".txt"
        synonymRejFileName = "Debug/pubchem_synonymRejections" + fileID_End_serialized + ".txt"
        outputFile = open(os.path.join(pyScript_directory , outputFileName),'w')
        outputErrors = open(os.path.join(pyScript_directory , debugFileName),'w')
        synonymRejectionFile = open(os.path.join(pyScript_directory , synonymRejFileName),'w')
        xmlFileInputPath = os.path.join(pyScript_directory , xmlFileInput)

        storedSyn = ""
        storedSyns = []
        outputSyns = []
        synStats = [0,0,0,0,0,0,0]
        searchSynBoolean = False
        synonymCountLimit = 50

        print("code executed")
        startTime = time.time()
        for event, elem in etree.iterparse(xmlFileInputPath):
                if elem.tag == "{http://www.ncbi.nlm.nih.gov}PC-CompoundType_id_cid": #print("switching to new CID: \t" + str(elem.text))
                        svalCount = 0
                        cid = elem.text.strip()
                if elem.tag == "{http://www.ncbi.nlm.nih.gov}PC-InfoData_value_sval": # reports everything from name to smiles
                        svalCount = svalCount + 1
                        if svalCount == 3:  #print("got iupac preferred name \t" + elem.text)
                                cidIUPAC = elem.text
                        if svalCount == 10: #print("got isomeric smiles \t" + elem.text)
                                cidISOSMILES = elem.text
                                searchSynBoolean = True
                if searchSynBoolean:
                        findSynCIDBoolean = True
                        synSearchBoolean = True
                        outputSyns = []
                        rejectSyns = []
                        rejectReasons = []
                        carryoverSyn = ""
                        while findSynCIDBoolean:
                                pcsLine = pubChemSyns.readline()
                                pcsLineSplit = pcsLine.split('\t')
                                cidSyn1 = pcsLineSplit[0].strip()
                                individualSyn = pcsLineSplit[1].strip() #print("searching for cid to match: \t" + str(cidSyn1) + "\t to real CID: \t" + str(cid))
                                if cidSyn1 == cid:
                                        findSynCIDBoolean = False #print("found first syn entry with matching cid: \t" + str(cidSyn1) + "\t to real CID: \t" + str(cid))
                                if (int(cidSyn1) > int(cid)): #print("synonym CID > current CID " + str(cidSyn1) + " " + str(cid) + ", so only one synonym was found for previous CID...")
                                        synSearchBoolean = False
                                        findSynCIDBoolean = False
                        storedSyns.append(individualSyn) #print("appended first synonym: \t" + individualSyn)
                        while synSearchBoolean:
                                pcsLine = pubChemSyns.readline()
                                pcsLineSplit = pcsLine.split('\t')
                                cidSyn2 = pcsLineSplit[0].strip()
                                individualSyn = pcsLineSplit[1].strip()
                                storedSyns.append(individualSyn)#print("compiling all synonyms with matching cid: \t" + str(cid) + "\t synonym: \t" + individualSyn)
                                # print("searching through matching cids: \t" + str(cidSyn2) + "\t to real CID: \t" + str(cid))
                                if cidSyn2 != cid:
                                        synSearchBoolean = False #print("reached next synonym set CID, popping last entry to next synonym set: \t" + individualSyn)                                
                        if len(storedSyns)>0:
                                #print("length greater than 0")
                                if len(storedSyns) < synonymCountLimit:#print("length less than limit")
                                        carryoverSyn = storedSyns.pop()
                                        for eachSyn in storedSyns:
                                                [synFilterBool,rejectionReasoning,rejectCode] = synFilter(eachSyn)
                                                if synFilterBool:
                                                        outputSyns.append(eachSyn)
                                                else:
                                                        rejectSyns.append(eachSyn)
                                                        rejectReasons.append(rejectCode)
                                        outputSyns.insert(0,cidIUPAC)
                                        output='|'.join(str(each) for each in outputSyns) + " \t root \t underivable \t " + cidISOSMILES + ",x" #print(output)
                                        outputFile.write(output+"\n")
                                        outputFile.flush()
                                        if len(outputSyns) < 6:
                                                synStats[len(outputSyns)] = synStats[len(outputSyns)] + 1
                                        else:
                                                synStats[6] = synStats[6] + 1
                                        outputSyns = []
                                else:
                                        outputErrors.write(cidIUPAC + "\t had too many usable, unfiltered synonyms: \t" + str(len(storedSyns)) + "\t limit was: \t" + str(synonymCountLimit)+"\n")
                                        outputErrors.flush() #print("too many usable syns: \t"+ str(len(storedSyns)))
                                        rejectSyns=rejectSyns + storedSyns
                                        rejectReasons.append("11: " + str(len(storedSyns)));
                        else:
                                outputErrors.write(cidIUPAC + "\t had no usable synonyms. amount was: \t" + str(len(storedSyns))+"\n")
                                outputErrors.flush() #print("error no usable syns stored from synonym bank: \t" + str(len(storedSyns)))
                                rejectSyns=rejectSyns + storedSyns
                                rejectReasons.append("12");
                        storedSyns = []
                        storedSyns.append(carryoverSyn) #print("carryover synonym moved: \t" + carryoverSyn)
                        searchSynBoolean = False
                        synonymRejectionReport = '|'.join(str(each) for each in rejectSyns) + "\t" + '|'.join(str(each) for each in rejectReasons) 
                        synonymRejectionFile.write(synonymRejectionReport + "\n") #print(synonymRejectionReport)
                        synonymRejectionFile.flush()
                elem.clear()    
        print("program terminated in: \t" + str(time.time()- startTime) + "\t with stats of syns: \t" + str(synStats))
        synStatsFile.write(str(fileID_Start) + "---" + str(fileID_End)+ "\t" + str(synStats) + "\n")
        synStatsFile.flush()
        multipleFileIterator = multipleFileIterator + 1
    
    
    
