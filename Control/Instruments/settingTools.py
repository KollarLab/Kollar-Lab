import re

def load_settings(inputfile):
    '''
    Helper function that reads in a file and turns it into a dictionary
    The file MUST have the following formatting standards:
    -[sections] that describe the main blocks/ modules of the instrument
    -[subsections] that describe subblocks such as channels
    -settings = value that match attributes of the section/subsection they are in
    Argurments:
        inputfile (str) : file where the settings are stored
    Returns:
        dictionary organized as follows:
        settings = {sectionA:{settings},sectionB:{subsectionB:{settings},etc.},etc.}
    '''

    #Regular expression that denote the style of the different setting labels
    sectionTag    = r'\[(\w*)\]'
    subsectionTag = r'\[\[(\w*)\]\]'
    settingTag    = r'(\S*)\s*=\s*(\S*)'

    confFile = open(inputfile,'r')
    line=confFile.readline()

    settings={}
    subset={}

    while True:
        #Determine the label associated with the line
        section    = re.match(sectionTag,line)
        subsection = re.match(subsectionTag,line)
        setting    = re.match(settingTag,line)

        if section:
            #Get the name of the section
            sectionName = section.group(1)
            #Read next line and check label
            line       = confFile.readline()
            section    = re.match(sectionTag,line)
            subsection = re.match(subsectionTag,line)
            setting    = re.match(settingTag,line)
            #Create dictionary for section entry
            settings[sectionName] = {}

            #Keep reading lines until we reach the next section
            while not section:
                #If we have a subsection, we need to go a layer deeper
                if subsection:
                    subset={}
                    #Get subsection name
                    subsectionName = subsection.group(1)
                    #Read next line
                    line    = confFile.readline()
                    setting = re.match(settingTag,line)
                    #Continue until we reach a non setting entry (section, subsection)
                    while setting:
                        subset[setting.group(1)] = setting.group(2)
                        #Read next line and check label
                        line       = confFile.readline()
                        setting    = re.match(settingTag,line)
                        section    = re.match(sectionTag,line)
                        subsection = re.match(subsectionTag,line)
                    #Store settings in the appropriate dict entry
                    settings[sectionName][subsectionName]=subset
                #If we just have settings, loop until we arrive to the next section
                elif setting:
                    while setting:
                        subset[setting.group(1)] = setting.group(2)
                        #Read next line and check label
                        line       = confFile.readline()
                        setting    = re.match(settingTag,line)
                        section    = re.match(sectionTag,line)
                        subsection = re.match(subsectionTag,line)
                    #Store settings in the section entry
                    settings[sectionName]=subset
                #If we reach the end of the file or something like that this should avoid infinite loops
                if not setting and not subsection:
                    break
        else:
            break

    return settings


def print_settings(settings):
    '''
    Helper function to display the settings stored in a dictionary
    Arguments:
        settings (str) : dictionary structured with nested dictionaries composed of sections, subsections and settings
    '''

    #Formatting the output
    sectionStyle    = '[{}]'
    subsectionStyle = '[[{}]]'
    settingStyle    = '{}={}'

    #Loop through all the entries in settins
    for key in settings.keys():
        #Print the section header
        print(sectionStyle.format(key))
        #Making sure that the section is actually populated
        if isinstance(settings[key],dict):
            for subkey in settings[key]:
                #Checks to see if we have a subsection
                if isinstance(settings[key][subkey],dict):
                    #Print subsection header
                    print(subsectionStyle.format(subkey))
                    #Loop over all the settings in the subsection
                    for setting in settings[key][subkey].keys():
                        print(settingStyle.format(setting,settings[key][subkey][setting]))
                else:
                    #Print the settings
                    print(settingStyle.format(subkey,settings[key][subkey]))

def save_settings(settings, savefile):
    '''
    Helper function to save settings to file, expects a dictionary with nested
    dictionaries coresponding to the various sections, subsections of the instrument settings
    Arguments:
        settings (dict) : dictionary containing all the relevant settings
        savefile (str) : path to file where settings should be saved
    '''

    #Formatting strings
    sectionStyle    = '[{}]\n'
    subsectionStyle = '[[{}]]\n'
    settingStyle    = '{}={}\n'

    saveTo=open(savefile,'w')
    #Loop through all the sections
    for key in settings.keys():
        saveTo.write(sectionStyle.format(key))
        #Check if section is populated
        if isinstance(settings[key],dict):
            for subkey in settings[key]:
                #Check if we have subsections
                if isinstance(settings[key][subkey],dict):
                    saveTo.write(subsectionStyle.format(subkey))
                    #Loop through all settings in subsection
                    for setting in settings[key][subkey].keys():
                        #Write settings to file
                        saveTo.write(settingStyle.format(setting,settings[key][subkey][setting]))
                else:
                    #Write settings to file
                    saveTo.write(settingStyle.format(subkey,settings[key][subkey]))

    saveTo.close()
            