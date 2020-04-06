import re

def load_settings(inputfile):
    sectionTag=r'\[(\w*)\]'
    subsectionTag=r'\[\[(\w*)\]\]'
    settingTag=r'(\S*)\s*=\s*(\S*)'

    confFile=open(inputfile,'r')
    settings={}
    subset={}
    line=confFile.readline()
    while True:
        section=re.match(sectionTag,line)
        subsection=re.match(subsectionTag,line)
        setting=re.match(settingTag,line)
        if section:
            sectionName=section.group(1)
            #print(sectionName)
            line=confFile.readline()
            section=re.match(sectionTag,line)
            subsection=re.match(subsectionTag,line)
            setting=re.match(settingTag,line)
            settings[sectionName]={}
            while not section:
                if subsection:
                    subset={}
                    subsectionName=subsection.group(1)
                    #print(subsectionName)
                    line=confFile.readline()
                    setting=re.match(settingTag,line)
                    while setting:
                        subset[setting.group(1)]=setting.group(2)
                        line=confFile.readline()
                        setting=re.match(settingTag,line)
                        section=re.match(sectionTag,line)
                        subsection=re.match(subsectionTag,line)
                    settings[sectionName][subsectionName]=subset
                elif setting:
                    while setting:
                        subset[setting.group(1)]=setting.group(2)
                        line=confFile.readline()
                        setting=re.match(settingTag,line)
                        section=re.match(sectionTag,line)
                        subsection=re.match(subsectionTag,line)
                    settings[sectionName]=subset
                if not setting and not subsection:
                    break
        else:
            break
    return settings


def print_settings(settings):
    sectionStyle='[{}]'
    subsectionStyle='[[{}]]'
    settingStyle='{}={}'

    for key in settings.keys():
        print(sectionStyle.format(key))
        if isinstance(settings[key],dict):
            for subkey in settings[key]:
                if isinstance(settings[key][subkey],dict):
                    print(subsectionStyle.format(subkey))
                    for setting in settings[key][subkey].keys():
                        print(settingStyle.format(setting,settings[key][subkey][setting]))
                else:
                    print(settingStyle.format(subkey,settings[key][subkey]))

def save_settings(settings, savefile):
    sectionStyle='[{}]\n'
    subsectionStyle='[[{}]]\n'
    settingStyle='{}={}\n'

    saveTo=open(savefile,'w')
    for key in settings.keys():
        saveTo.write(sectionStyle.format(key))
        if isinstance(settings[key],dict):
            for subkey in settings[key]:
                if isinstance(settings[key][subkey],dict):
                    saveTo.write(subsectionStyle.format(subkey))
                    for setting in settings[key][subkey].keys():
                        saveTo.write(settingStyle.format(setting,settings[key][subkey][setting]))
                else:
                    saveTo.write(settingStyle.format(subkey,settings[key][subkey]))
    saveTo.close()
            