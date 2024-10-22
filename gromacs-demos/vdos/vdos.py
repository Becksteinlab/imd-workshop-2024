# %%
import ctypes as ct
import MDAnalysis as mda
from scipy.fft import fft
import numpy as np

# %%
# ctype data structures for rotatable bonds
class t_rotBond(ct.Structure):
    '''datatype describing a rotatable bond and its satellites'''
    _fields_ = (("bondAtomIndices", ct.c_int32*2),
                ("sat0AtomIndices", ct.c_int32*4),
                ("nSat0", ct.c_int32),
                ("sat3AtomIndices", ct.c_int32*4),
                ("nSat3", ct.c_int32),
                ("sumInertia", ct.c_double),
                ("sumLogInertia", ct.c_double),
                ("wOmegaBuffer", ct.POINTER(ct.c_float))
                )

class t_residue(ct.Structure):
    '''datatype with data and pointers for residue/molecule)'''
    _fields_ = (# constant data
                ("nAtoms", ct.c_int32),
                ("resMass", ct.c_float),
                ("nRotBonds", ct.c_int32),
                ("atomMasses", ct.POINTER(ct.c_float)),
                ("atomMassesX3", ct.POINTER(ct.c_float)),
                ("rotBonds", ct.POINTER(t_rotBond)),
                # transient data
                ("atomCrdMol", ct.POINTER(ct.c_float)),
                ("atomVelMol", ct.POINTER(ct.c_float)),
                ("inertia", ct.c_float*3),
                ("prevAxes", ct.c_float*9),
                ("angMomLab", ct.c_float*3),
                ("angMomMol", ct.c_float*3),
                ("omegaMol", ct.c_float*3),
                ("omegaLab", ct.c_float*3),
                ("COMposBuffer", ct.POINTER(ct.c_float)),
                ("COMposCurrent", ct.POINTER(ct.c_float)),
                ("COMposCorrRef", ct.POINTER(ct.c_float)),
                ("COMvelBuffer", ct.POINTER(ct.c_float)),
                ("COMvelCurrent", ct.POINTER(ct.c_float)),
                ("COMvelCorrRef", ct.POINTER(ct.c_float)),
                ("wOmegaBuffer", ct.POINTER(ct.c_float)),
                ("wOmegaCurrent", ct.POINTER(ct.c_float)),
                ("wOmegaCorrRef", ct.POINTER(ct.c_float)),
                # accumulating data
                ("sumInertia", ct.c_double*3),
                ("sumLogInertia", ct.c_double*3),
                ("totCorr",ct.POINTER(ct.c_double)),
                ("trCorr",ct.POINTER(ct.c_double)),
                ("rotCorr",ct.POINTER(ct.c_double)),
                ("rotBondCorr",ct.POINTER(ct.c_double)),
                ("propCnt", ct.c_int32),
                ("corrCnt", ct.c_int32),
                # offset in arrays above
                ("offset", ct.c_int32)
                )

class t_residueList(ct.Structure):
    '''datatype describing a list of residues (e.g., for each residue in a selection)'''
    _fields_ = (("residues", ct.POINTER(t_residue)),
                ("nResidues", ct.c_int32),
                ("totCorr",ct.POINTER(ct.c_double)),
                ("trCorr",ct.POINTER(ct.c_double)),
                ("rotCorr",ct.POINTER(ct.c_double)),
                ("rotBondCorr",ct.POINTER(ct.c_double)),
                ("inertia", ct.c_double*3),
                ("logInertia", ct.c_double*3),
                ("rotBondInertia", ct.c_double),
                ("logRotBondInertia", ct.c_double),
                ("postProcessed", ct.c_int32)
                )
    
class t_MDinfo(ct.Structure):
    '''datatype with data and pointers for MD simulation'''
    _fields_ = (# number of floats in per frame in array
                ("frameSize", ct.c_int32),
                # single set of atomic coordinates for selection
                ("atomCrd", ct.POINTER(ct.c_float)),
                # number of floats in atomVelBuffer
                ("bufferSize", ct.c_int32),
                # nCorr sets of atomic velocities for selection
                ("atomVelBuffer", ct.POINTER(ct.c_float)),
                # pointer to current set of atomic velocities in buffer
                ("atomVelCurrent", ct.POINTER(ct.c_float)),
                # pointer to first set of atomic velocities in current buffer
                ("atomVelCorrRef", ct.POINTER(ct.c_float)),
                # nCorr sets of atomic velocities for selection
                ("nCorr", ct.c_int32),
                # number of atoms in selection
                ("nAtomsSel", ct.c_int32)
                )

# %%
#load the shared library with C routines
vdosLib = ct.cdll.LoadLibrary("libvdos.so")

# %%
# allocate residues in residueList
vdosLib.allocResidueList.argtypes = [
    # list of residues in selection
    ct.POINTER(t_residueList),
    # number of residues in selection
    ct.c_int32,
    # number of correlation times to allocate arrays
    ct.c_int32
]
vdosLib.allocResidueList.restype = ct.c_int32

# allocate MDinfo
vdosLib.allocMDinfo.argtypes = [
    # list of residues in selection
    ct.POINTER(t_MDinfo),
    # number of correlation times
    ct.c_int32,
    # number of atoms in selection
    ct.c_int32
]
vdosLib.allocMDinfo.restype = ct.c_int32

# allocate all arrays needed for residue and set constant pointers
vdosLib.allocResidue.argtypes = [
    # pointer to residue
    ct.POINTER(t_residue),
    # number of atoms in residue
    ct.c_int32,
    # # list of atomic masses for residue
    np.ctypeslib.ndpointer(dtype=np.float32, ndim=1, flags='C_CONTIGUOUS'),
    # residue mass
    ct.c_float,
    # number of atoms in selection
    # (needed for offsets in atomic velocity buffer)
    ct.c_int32,
    # number of correlation times to allocate arrays
    ct.c_int32
]
vdosLib.allocResidue.restype = ct.c_int32

# set pointers to atom coordinates and velocities for each residue
vdosLib.setArrayIndexOffsets.argtypes = [
    # pointer to residue list
    ct.POINTER(t_residueList),
    # pointer to MDinfo
    ct.POINTER(t_MDinfo)
]
vdosLib.setArrayIndexOffsets.restype = ct.c_int32

#find rotatable bonds in all residues
vdosLib.getRotBonds.argtypes = [
    # sets of residues
    ct.POINTER(t_residueList),
    # number of residues
    ct.c_int32,
    # list of atom indices for dihedrals
    np.ctypeslib.ndpointer(dtype=np.int32, ndim=2, flags='C_CONTIGUOUS'),
    # number of dihedrals
    ct.c_int32,
    # pair of atom indices indicating the first and last atom index of each residue
    np.ctypeslib.ndpointer(dtype=np.int32, ndim=2, flags='C_CONTIGUOUS'),
    # number of correlation times to allocate arrays
    ct.c_int32
]
vdosLib.getRotBonds.restype = ct.c_int32

#process simulation time step
vdosLib.processStep.argtypes = [
    #time step
    ct.c_int32,
    #structure to store MD atomic corrdinates and buffer of velocities
    ct.POINTER(t_MDinfo),
    #atomic coordinates for selection
    np.ctypeslib.ndpointer(dtype=np.float32, ndim=2, flags='C_CONTIGUOUS'),
    #atomic velocities for selection
    np.ctypeslib.ndpointer(dtype=np.float32, ndim=2, flags='C_CONTIGUOUS'),
    #residue list to store residue properties
    ct.POINTER(t_residueList),
    #number of correlation times to allocate arrays
    ct.c_int32
]
vdosLib.processStep.restype = ct.c_int32

vdosLib.copyResidueListData.argtypes = [
    #pointer to residue list copy (dst)
    ct.POINTER(t_residueList),
    #pointer to residue list (src))
    ct.POINTER(t_residueList),
    #number of correlation times
    ct.c_int32
]
vdosLib.copyResidueListData.restype = ct.c_int32

# normalize data in single residue
vdosLib.normResData.argtypes = [
    ct.POINTER(t_residue),
    #number of correlation times
    ct.c_int32
]
vdosLib.normResData.restype = ct.c_int32

# normalize data for all residues and sum up for residue list
vdosLib.sumData.argtypes = [
    #time step
    ct.POINTER(t_residueList),
    #number of correlation times
    ct.c_int32
]
vdosLib.sumData.restype = ct.c_int32

# %%
class vdos:
    def __init__(self,sel,nCorr):
        self.sel = sel
        self.nCorr = nCorr
        self.nRes = sel.residues.n_residues
        # time and frequency axes for correlation functions and VDoS
        self.tau = np.zeros(nCorr, dtype = np.float32)
        self.wavenumber = np.zeros(nCorr, dtype = np.float32)
        # numpy arrays for average (index 0) and per residue VACF & VDoS (in 3D-2PT for each voxel)
        # -> total VACF / VDoS (1D)
        # -> translation VACF / VDoS (3D)
        # -> rotation VACF / VDoS (3D)
        # -> rotatable bonds VACF / VDoS (1D)
        self.totVACF     = np.zeros((self.nRes+1,    nCorr), dtype = np.float32)
        self.totVDoS     = np.zeros((self.nRes+1,    nCorr), dtype = np.float32)
        self.trVACF      = np.zeros((self.nRes+1, 3, nCorr), dtype = np.float32)
        self.trVDoS      = np.zeros((self.nRes+1, 3, nCorr), dtype = np.float32)
        self.rotVACF     = np.zeros((self.nRes+1, 3, nCorr), dtype = np.float32)
        self.rotVDoS     = np.zeros((self.nRes+1, 3, nCorr), dtype = np.float32)
        self.rotBondVACF = np.zeros((self.nRes+1,    nCorr), dtype = np.float32)
        self.rotBondVDoS = np.zeros((self.nRes+1,    nCorr), dtype = np.float32)
        # initialize residueList and MDinfo
        self.residueList, self.MDinfo, self.residueListCopy = self.prep()
        self.residueList.postProcessed = 0
        self.residueListCopy.postProcessed = 0

    def prep(self):
        '''construct list of residues for selection'''
        residueList = t_residueList()
        error = vdosLib.allocResidueList(
            ct.pointer(residueList),
            ct.c_int(self.nRes),
            ct.c_int(self.nCorr)
        )
        MDinfo = t_MDinfo()
        error = vdosLib.allocMDinfo(
            ct.pointer(MDinfo),
            ct.c_int(self.nCorr),
            ct.c_int(self.sel.n_atoms)
        )
        r = 0
        for res in self.sel.residues:
            error = vdosLib.allocResidue(
                ct.pointer(residueList.residues[r]),
                ct.c_int(len(res.atoms)),
                res.atoms.masses.astype(np.float32),
                ct.c_float(res.mass.astype('float32')),
                ct.c_int(self.sel.n_atoms),
                ct.c_int(self.nCorr)
            )
            if error != 0:
                print(f'ERROR reported by \'vdosLib.allocResidue\'\n')
            r += 1
        vdosLib.setArrayIndexOffsets(ct.pointer(residueList), ct.pointer(MDinfo))
        dihed = self.sel.intra_dihedrals.indices.astype(np.int32)
        resAtomRangeList = np.zeros((len(self.sel.residues),2), dtype = np.int32)
        i = 0
        for res in self.sel.residues:
            resAtomRangeList[i][0] = np.amin(res.atoms.indices)
            resAtomRangeList[i][1] = np.amax(res.atoms.indices)
            i += 1
        error = vdosLib.getRotBonds(
            ct.pointer(residueList),
            ct.c_int(self.nRes),
            dihed,
            ct.c_int(len(dihed)),
            resAtomRangeList,
            self.nCorr
        )
        if error != 0:
            print(f'ERROR reported by \'vdosLib.getRotBonds\'\n')

        # the copy is meant for intermediate postprocessing, i.e., live monitoring
        residueListCopy = t_residueList()
        error = vdosLib.allocResidueList(
            ct.pointer(residueListCopy),
            ct.c_int(self.nRes),
            ct.c_int(self.nCorr)
        )
        r = 0
        for res in self.sel.residues:
            error = vdosLib.allocResidue(
                ct.pointer(residueListCopy.residues[r]),
                ct.c_int(len(res.atoms)),
                res.atoms.masses.astype(np.float32),
                ct.c_float(res.mass.astype('float32')),
                ct.c_int(self.sel.n_atoms),
                ct.c_int(self.nCorr)
            )
            if error != 0:
                print(f'ERROR reported by \'vdosLib.allocResidue\'\n')
            r += 1
        vdosLib.setArrayIndexOffsets(ct.pointer(residueListCopy), ct.pointer(MDinfo))
        error = vdosLib.getRotBonds(
            ct.pointer(residueListCopy),
            ct.c_int(self.nRes),
            dihed,
            ct.c_int(len(dihed)),
            resAtomRangeList,
            self.nCorr
        )
        if error != 0:
            print(f'ERROR reported by \'vdosLib.getRotBonds\'\n')
        return residueList, MDinfo, residueListCopy
    
    def single_frame(self,tStep,time):
        if tStep < self.nCorr:
            self.tau[tStep] = time
        error=vdosLib.processStep(
            ct.c_int(tStep),
            ct.pointer(self.MDinfo),
            self.sel.atoms.positions.astype(np.float32),
            self.sel.atoms.velocities.astype(np.float32),
            ct.pointer(self.residueList),
            ct.c_int(self.nCorr)
        )

    def copyResidueList(self):
        '''copy residueList to residueListCopy'''
        vdosLib.copyResidueListData(
            ct.pointer(self.residueListCopy), 
            ct.pointer(self.residueList), 
            ct.c_int(self.nCorr)
        )
        self.residueListCopy.postProcessed = 0        

    def postProcess(self,residueList, mode = "all", resIdx = 0):
        '''post-process correlation functions'''
        if residueList.postProcessed == 0:
            period = (self.tau[1] - self.tau[0]) * (2 * self.nCorr - 1)
            wn0 = (1.0 / period) * 33.35641
            self.wavenumber = np.arange(0,self.nCorr) * wn0
            if mode =="single":
                vdosLib.normResData(ct.pointer(residueList.residues[resIdx]), ct.c_int(self.nCorr))
            elif mode == "all" or mode == "total" or mode == "total+single":
                vdosLib.sumData(ct.pointer(residueList), ct.c_int(self.nCorr))
                self.totVACF[0] = residueList.totCorr[0:self.nCorr]
                self.symFT(self.totVACF[0], self.totVDoS[0])
                for i in range(3):
                    self.trVACF[0][i] = residueList.trCorr[i:3*self.nCorr:3]
                    self.symFT(self.trVACF[0][i], self.trVDoS[0][i])
                    self.rotVACF[0][i] = residueList.rotCorr[i:3*self.nCorr:3]
                    self.symFT(self.rotVACF[0][i], self.rotVDoS[0][i])
                self.rotBondVACF[0] = residueList.rotBondCorr[0:self.nCorr]
                self.symFT(self.rotBondVACF[0], self.rotBondVDoS[0])
            else:
                print('ERROR: unknown mode in postProcess')
                raise ValueError
            if mode == "all":
                for i in range(1,self.nRes+1):
                    self.totVACF[i] = residueList.residues[i-1].totCorr[0:self.nCorr]
                    self.symFT(self.totVACF[i], self.totVDoS[i])
                    for j in range(3):
                        self.trVACF[i][j] = residueList.residues[i-1].trCorr[j:3*self.nCorr:3]
                        self.symFT(self.trVACF[i][j], self.trVDoS[i][j])
                        self.rotVACF[i][j] = residueList.residues[i-1].rotCorr[j:3*self.nCorr:3]
                        self.symFT(self.rotVACF[i][j], self.rotVDoS[i][j])
                    self.rotBondVACF[i] = residueList.residues[i-1].rotBondCorr[0:self.nCorr]
                    self.symFT(self.rotBondVACF[i], self.rotBondVDoS[i])
            elif mode == "single" or mode == "total+single":
                if mode == "single":
                    i = 0
                else:
                    i = 1
                self.totVACF[i] = residueList.residues[resIdx].totCorr[0:self.nCorr]
                self.symFT(self.totVACF[i], self.totVDoS[i])
                for j in range(3):
                    self.trVACF[i][j] = residueList.residues[resIdx].trCorr[j:3*self.nCorr:3]
                    self.symFT(self.trVACF[i][j], self.trVDoS[i][j])
                    self.rotVACF[i][j] = residueList.residues[resIdx].rotCorr[j:3*self.nCorr:3]
                    self.symFT(self.rotVACF[i][j], self.rotVDoS[i][j])
                self.rotBondVACF[i] = residueList.residues[resIdx].rotBondCorr[0:self.nCorr]
                self.symFT(self.rotBondVACF[i], self.rotBondVDoS[i])
            residueList.postProcessed = 1

    def symFT(self,input,output):
        '''symmetrize, Fourier transform and return first half of spectrum'''
        data = np.array(input, dtype = np.float32)
        nData=len(data)
        tmp = np.zeros(2 * nData - 1, dtype = np.float32)
        for i in range(nData):
            tmp[i] = data[i]
        for i in range(1,nData):
            j = 2 * nData - i - 1
            tmp[j] = tmp[i]
        tmp = fft(tmp)
        for i in range(nData):
            output[i] = tmp[i].real
    
    def outputGeometry(self,outputFileName,residueList):
        totMass = 0.0
        for i in range(self.nRes):
            totMass += residueList.residues[i].resMass
        averMass = totMass / self.nRes
        averInertia = residueList.inertia[0:3]
        averLogInertia = residueList.logInertia[0:3]
        averRotBondCount = 0
        for i in range(self.nRes):
            averRotBondCount += residueList.residues[i].nRotBonds
        averRotBondCount /= self.nRes
        if averRotBondCount > 0:
            averRotBondInertia = residueList.rotBondInertia / averRotBondCount
            averLogRotBondInertia = residueList.logRotBondInertia / averRotBondCount
        else:
            averRotBondInertia = 0.0
            averLogRotBondInertia = 0.0
        outFile = open(outputFileName,"w")
        outFile.write("#Average Mass (g/mol):\n%20.6e\n" % averMass)
        outFile.write("#Average Inertia (g/mol*A^2):\n%20.6e %20.6e %20.6e\n" % (averInertia[0],averInertia[1],averInertia[2]))
        outFile.write("#Average Log(Inertia):\n%20.6e %20.6e %20.6e\n" % (averLogInertia[0],averLogInertia[1],averLogInertia[2]))
        outFile.write("#Average Rotatable Bond Inertia (g/mol*A^2): (%d)\n%20.6e\n" % (averRotBondCount,averRotBondInertia))
        outFile.write("#Average Rotatable Bond Log(Inertia): (%d)\n%20.6e\n" % (averRotBondCount,averLogRotBondInertia))
        outFile.close()
    
    def outputVACF(self,outputFileName):
        outFile = open(outputFileName,"w")
        outFile.write("%-20s %-20s %-20s %-20s %-20s %-20s %-20s %-20s %-20s\n" % ("#Time (ps)","Translation x","Translation y","Translation z","Rotation x","Rotation y","Rotation z","Rotatable Bonds","Total"))
        for i in range(self.nCorr):
            outFile.write("%-20.6e %-20.6e %-20.6e %-20.6e %-20.6e %-20.6e %-20.6e %-20.6e %-20.6e\n" % (self.tau[i],self.trVACF[0][0][i],self.trVACF[0][1][i],self.trVACF[0][2][i],self.rotVACF[0][0][i],self.rotVACF[0][1][i],self.rotVACF[0][2][i],self.rotBondVACF[0][i],self.totVACF[0][i]))
        outFile.close()
        outFile.close()

    def outputVDoS(self,outputFileName):
        outFile = open(outputFileName,"w")
        outFile.write("%-20s %-20s %-20s %-20s %-20s %-20s %-20s %-20s %-20s\n" % ("#Wavenumber (cm^-1)","Translation x","Translation y","Translation z","Rotation x","Rotation y","Rotation z","Rotatable Bonds","Total"))
        for i in range(self.nCorr):
            outFile.write("%-20.6e %-20.6e %-20.6e %-20.6e %-20.6e %-20.6e %-20.6e %-20.6e %-20.6e\n" % (self.wavenumber[i],self.trVDoS[0][0][i],self.trVDoS[0][1][i],self.trVDoS[0][2][i],self.rotVDoS[0][0][i],self.rotVDoS[0][1][i],self.rotVDoS[0][2][i],self.rotBondVDoS[0][i],self.totVDoS[0][i]))
        outFile.close()