class ParseToJSON:
	'''
	This function parses input data into json format
	'''
	def __init__(self,list_of_files):
		'''
		This function reads the input data into a dictionary, converts it to json, and prints it to an outfile
		'''
		outfile = open('23andMe.json', 'w')
		self.Data = {}
		for i in list_of_files:
				self.Data[i] = self.readInFile(i)
		self.Data = self.convertToJSON(self.Data)
		outfile.write(str(self.Data))
		outfile.close()

	def readInFile(self,infile):
		'''
		This function reads in the data from a file into a dictionary
		'''
		f = open(infile, 'r')
		rsid = {}
		for i in f.readlines():
			if '#' not in i:
				line = i.strip().split('\t')
				rsid[line[0]] = line[1:]
		f.close()
		return rsid

	def convertToJSON(self,data):
		'''
		This function converts the data structure (dictionary) into a json string
		'''
		out = JSONEncoder().encode(data)
		return out
