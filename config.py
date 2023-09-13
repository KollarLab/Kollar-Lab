import site
import os
import sys
from distutils.sysconfig import get_python_lib

#site_directory = site.getusersitepackages()
#
#site_directory = r'C:\Users\Kollarlab\Anaconda3\Lib\site-packages'

site_directory = get_python_lib()

package_directory = os.getcwd()

filepath = os.path.join(site_directory,"ControlCode.pth")

print(site_directory)
print(package_directory)
print(filepath)

f = open(filepath,"w")
f.write(package_directory)