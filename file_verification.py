import os,sys
    
def file_verification(path,filename,directory):
	
	if type(directory)==int or type(directory)==float:
		directory = str(directory)
	elif type(directory)==str:
		pass	
	else:
		raise(sys.exc_info()[0])
		
	if os.path.isdir(path):
		path = os.path.join(path,directory) 
		if os.path.isdir(path):
			path = os.path.join(path,filename)
			if os.path.isfile(path):
				os.remove(path)
		else:
			os.mkdir(path)
	else:
		os.mkdir(path)
		os.mkdir(os.path.join(path,directory))
