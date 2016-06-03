bashCommand = "nohup nice matlab -nodisplay -nosplash -r main > test_output.txt /dev/null 2>&1 &"
import subprocess
process = subprocess.Popen(bashCommand.split(), stdout=subprocess.PIPE)
output = process.communicate()[0]
