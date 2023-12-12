Spyder 5:

1. The safest way to launch Spyder in a different virtual environment is from Anaconda Navigator. You might see different Spyder apps that have the name of the virtual environment you want to launch, this is a trap, they won't necessarily launch in that environment. You must either use Anaconda Navigator or activate the environment from the command line to be sure it's opened correctly.

2. **IMPORTANT** When a new instance of Spyder 5 is launched, the default behavior is for individual files to store variables in individual namespaces. If a variable is defined in one file, the other files won't be able to access it, which breaks our measurement scripts. To turn it back on check the box in the pathway below: 
Tools -> Preferences -> Run in console namespace instead of an empty one

3. How to change Spyder's default to make figures in a separate window: 
Tools -> Preferences -> iPython Console -> Graphics -> Change backend from Inline to Automatic

4. How to turn on the shift + enter shortcut for running highlighted lines of code:
Tools -> Preferences -> Keyboard Shortcuts -> Scroll down to run selection and change default option (usually F9) to shift+enter
