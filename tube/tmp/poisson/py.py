import subprocess
import os
import tempfile
import logfile
processes = []
for file in ['lambda2_0_t=0.0336825_tracer.png','lambda2_0_t=3.2526_noLambda2.png',    'lambda2_0_t=3.2526ux.png']:
    f = tempfile.TemporaryFile()
    p = subprocess.Popen(['ls',file],stdout=f)
    processes.append((p, f))

for p, f in processes:
    p.wait()
    f.seek(0)
    logfile.write(f.read())
    f.close()
