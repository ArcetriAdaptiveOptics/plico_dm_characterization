'''
Created on 1 lug 2024

@author: cselmi
'''
import os
from plico_dm_characterization.configuration import config

def main():
    if os.path.exists(config.ROOT_FOLDER) == False:
        os.makedirs(config.ROOT_FOLDER)
    
    if os.path.exists(config.MODALAMPLITUDE_ROOT_FOLDER) == False:
        os.makedirs(config.MODALAMPLITUDE_ROOT_FOLDER)
    if os.path.exists(config.MODALBASE_ROOT_FOLDER) == False:
        os.makedirs(config.MODALBASE_ROOT_FOLDER)
    if os.path.exists(config.COMMANDHISTORY_ROOT_FOLDER) == False:
        os.makedirs(config.COMMANDHISTORY_ROOT_FOLDER)
    if os.path.exists(config.IFFUNCTIONS_ROOT_FOLDER) == False:
        os.makedirs(config.IFFUNCTIONS_ROOT_FOLDER)
    if os.path.exists(config.FLAT_ROOT_FOLD) == False:
        os.makedirs(config.FLAT_ROOT_FOLD)
    if os.path.exists(config.OPD_SERIES_ROOT_FOLD) == False:
        os.makedirs(config.OPD_SERIES_ROOT_FOLD)
    if os.path.exists(config.LINEARITY_ROOT_FOLD) == False:
        os.makedirs(config.LINEARITY_ROOT_FOLD)

    print('All yours folder have been created here: %s' %config.ROOT_FOLDER)

if __name__ == "__main__":
    main()