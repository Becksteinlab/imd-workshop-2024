# %%
from imdclient.IMDREADER import IMDReader
import MDAnalysis as mda
import align as al
import logging

# %%
logger = logging.getLogger("imdclient.IMDClient")
file_handler = logging.FileHandler("imdreader.log")
formatter = logging.Formatter(
    "%(asctime)s - %(name)s - %(levelname)s - %(message)s"
)
file_handler.setFormatter(formatter)
logger.addHandler(file_handler)
logger.setLevel(logging.DEBUG)

# %%
u = mda.Universe("mda.tpr", "imd://localhost:8888", buffer_size = 100*1024*1024)
protein = u.select_atoms("protein")
align = al.align(u,protein)

i = 0
with mda.Writer("imd/protein-pbc-aligned.trr", protein.n_atoms) as w:
    for ts in u.trajectory:
        align.single_frame()
        if i % 10 == 0:
            w.write(protein)
        i += 1

# %%
logger.info(f"Parsed {i} frames")


