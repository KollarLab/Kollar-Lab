import clr
clr.AddReference('pydll')

import ClassLib

testobj=ClassLib.class1()

print(testobj.Function1())
testobj.Function2()