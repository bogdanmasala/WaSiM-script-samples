'''
Code for deleting unneeded files from WaSiM output folders.

'''

#import libraries
import os

#get and set path of working directory
path = os.getcwd() 
os.chdir(working directory)

#set directory paths
output_directory = 'output directory' #where output folders are located

##############################################################################################################################################################################
def delete_files(output_directory):
    
    output_folder_names = []
    i=0
    for i in range(Nstart, Nend):
        output_folder_names.append('output folder base name' + str(i+1))

    for a in range(Nstart, Nend):
        output_folder = os.path.join(output_directory, output_folder_names[a])
        #lines below if deleting output files after each model run
        for delete_file in os.listdir(output_folder):
            if delete_file.startswith('file start') and delete_file.endswith('file ending'):
            #if not delete_file.startswith(('25_Martell_','glc_','glid_','Snow_','special_','ssto_','sstodem_','table_')):
                os.unlink(os.path.join(output_folder, delete_file))
                
    #end function
##############################################################################################################################################################################

#run function
delete_files(output_directory)    
