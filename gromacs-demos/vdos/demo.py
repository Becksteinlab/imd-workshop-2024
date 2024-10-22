from imdclient.IMDREADER import IMDReader
import numpy as np
import MDAnalysis as mda
import logging
import dynplot as dyn
import vdos as vd
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

# %%
u = mda.Universe("mda.tpr", "imd://localhost:8888",buffersize = 10*1024*1024)
sel = u.select_atoms("resname SOL")
vdos = vd.vdos(sel,200)

vd.vdosLib.omp_set_num_threads(2)
tStep = 0
for ts in u.trajectory:
    vdos.single_frame(tStep,ts.time)
    if tStep == 0:
        d = dyn.dynamicPlot([0,1000],"wavenumber (cm^-1)",[0,12000],"VDoS (arb. u.)","")
    elif tStep % 200 == 0:
        vdos.copyResidueList()
        vdos.postProcess(vdos.residueListCopy,mode = "total+single")
        total = np.array(vdos.totVDoS[0])
        single = np.array(vdos.totVDoS[1])
        d.update(vdos.wavenumber, total, single)
    tStep += 1
# %%
logger.info(f"Parsed {tStep} frames")


