# %%
from imdclient.IMDREADER import IMDReader
import numpy as np
import MDAnalysis as mda
import logging
import dynplot as dyn
import time as t

# %%
logger = logging.getLogger("imdclient.IMDClient")
file_handler = logging.FileHandler("imdreader.log")
formatter = logging.Formatter(
    "%(asctime)s - %(name)s - %(levelname)s - %(message)s"
)
file_handler.setFormatter(formatter)
logger.addHandler(file_handler)
logger.setLevel(logging.INFO)

timeWindow = 5.0
dtIMD = 0.010
nTimesPlot = (int)(timeWindow/dtIMD)

# %%
u = mda.Universe("mda.tpr", "imd://localhost:8888", buffer_size = 100*1024*1024)
Nter = u.select_atoms("resid 1 and name CA")
Cter = u.select_atoms("resid 129 and name CA")
GLU35 = u.select_atoms("resid 35 and name CD")
ASP52 = u.select_atoms("resid 52 and name CG")
time = np.zeros(nTimesPlot, dtype=np.float32)
dist1 = np.zeros(nTimesPlot, dtype=np.float32)
dist2 = np.zeros(nTimesPlot, dtype=np.float32)

i = 0
for ts in u.trajectory:
    idx = i % nTimesPlot
    time[idx] = ts.time
    dist1[idx] = np.linalg.norm(Nter.atoms[0].position - Cter.atoms[0].position)
    dist2[idx] = np.linalg.norm(GLU35.atoms[0].position - ASP52.atoms[0].position)
    if i == 0:
        d = dyn.dynamicPlot([0,timeWindow],"time (ps)",[18,22],[3,10],"distance (A)",["Nter-Cter","Glu35-Asp52"])
    elif i % nTimesPlot ==0:
        d.update2(time, dist1, dist2, [time[0], time[0]+timeWindow])
    elif i % 1 == 0:
        d.update(time[:idx], dist1[:idx], dist2[:idx])
    i += 1
# %%
logger.info(f"Parsed {i} frames")


