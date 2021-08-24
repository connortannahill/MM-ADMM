import time
import os
import subprocess
bashCmd="./mesh.exe"
start = time.perf_counter()
for i in range(10):
    #process = subprocess.Popen(bashCmd, stdout=subprocess.PIPE)
    os.system(bashCmd)
end = time.perf_counter()
print(f"Total time: {end - start:0.4f} seconds.")
print(f"Average time: {(end-start)/10:0.4f} seconds.")