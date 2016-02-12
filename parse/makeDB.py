import os, sys
import sqlite3
import json
from json import JSONEncoder

class ParseToDB:
	'''
	This class loads the input files into a sqlite database
	'''
	def __init__(self,list_of_files):
		'''
		This function parses the 23andMe raw data file to a sqlite database
		'''
		dbName = '<yourDB_name>.db'
		try:
			self.createDB(dbName,list_of_files)
		except sqlite3.OperationalError:
			print "Database", dbName, "already exists. Moving on..."
			return None
		for i in list_of_files:
				self.readInFilesToDB(i,dbName)
		
	def createDB(self,name,list_of_files):
		'''
		This function creates a database
		'''
		DB = sqlite3.connect(name)
		cursor = DB.cursor()
		for i in list_of_files:
			TableName = str(i.split('.')[0])
			cursor.execute('create table '+TableName+' (rsid TEXT, chr TEXT, pos TEXT, geno TEXT)')
		DB.close()
	
	def readInFilesToDB(self,infile,dbName):
		'''
		This function reads in the 23andMe raw data file into a dictionary
		'''
		f = open(infile,'r')
		DB = sqlite3.connect(dbName)
		cursor = DB.cursor()
		for i in f.readlines():
			if '#' not in i:
				line = i.strip().split('\t')
				TableName = str(infile.split('.')[0])
				cursor.execute("insert into "+TableName+" values(?,?,?,?)", (line[0],line[1:][0],line[1:][1],line[1:][2]))
		DB.commit()
		DB.close()
		f.close()


