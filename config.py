import site
import os
import sys

site_directory = site.getusersitepackages()

site_directory = r'C:\Users\Martin\AppData\Local\Programs\Python\Python39\Lib\site-packages'
package_directory = os.getcwd()

filepath = os.path.join(site_directory+"\kollar_measurement.pth")

print(site_directory)
print(package_directory)
print(filepath)

f = open(filepath,"w")
f.write(package_directory)
f.write('\n')