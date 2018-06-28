# Get the name of files in a directory and its subdirectories
# Author: Wei Tian <Wei.Tian@Schneider-Electric.com>
# Last Update: Thu, 28 June 2018
'''
A FUNCTION TO READ THE FIRTRAN FILE NAMES IN DIRECTORY
AND WRITE IT TO FILE NAMED FILNAM
'''
def getFilesNames(directory):
    from os import walk
    f = []
    for (dirpath, dirnames, filenames) in walk(directory):
        for filename in filenames:
            print (filename)
            if ('.f90' in filename):
                f.append(filename)
            else:
                continue
    return (f)

def writeSourceFile(filenames,writefilnam):
    with open (writefilnam,'w') as f:
        f.write('# Get the name of files in a directory and its subdirectories\n')
        f.write('# Author: Wei Tian <Wei.Tian@Schneider-Electric.com>\n')
        f.write('FILES=')
        for index, filename in enumerate(filenames):
            if((index+1)%5 == 0):
                f.write('%s    \ \n\t' %filename)
            else:
                f.write('%s ' %filename)
        

if __name__ == '__main__':
    readdirectory = '../src/isat'
    writefilnam = '../compile/isat/SOURCES.mk'
    filenames = getFilesNames(readdirectory)
    writeSourceFile(filenames,writefilnam)
