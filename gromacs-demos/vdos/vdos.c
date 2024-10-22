#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_eigen.h>

typedef float *vector;

void vecCopy(vector a, vector res) {
    res[0]=a[0];
    res[1]=a[1];
    res[2]=a[2];
}

void vecSub(vector a, vector b, vector res) {
    res[0]=a[0]-b[0];
    res[1]=a[1]-b[1];
    res[2]=a[2]-b[2];
}

void vecSubInPlace(vector a, vector b) {
    a[0]-=b[0];
    a[1]-=b[1];
    a[2]-=b[2];
}

void vecAdd(vector a, vector b, vector res) {
    res[0]=a[0]+b[0];
    res[1]=a[1]+b[1];
    res[2]=a[2]+b[2];
}

void vecAddInPlace(vector a, vector b) {
    a[0]+=b[0];
    a[1]+=b[1];
    a[2]+=b[2];
}

float vecDot(vector a, vector b) {
    return a[0]*b[0]+a[1]*b[1]+a[2]*b[2];
}

float vecNorm(vector a) {
    return sqrt(vecDot(a,a));
}

void vecNormalize(vector a, vector res) {
    float norm;
    norm=vecNorm(a);
    res[0]=a[0]/norm;
    res[1]=a[1]/norm;
    res[2]=a[2]/norm;
}

void vecNormalizeInPlace(vector a) {
    float norm;
    norm=vecNorm(a);
    a[0]/=norm;
    a[1]/=norm;
    a[2]/=norm;
}

void vecScale(float s, vector a, vector res) {
    res[0]=a[0]*s;
    res[1]=a[1]*s;
    res[2]=a[2]*s;
}

void vecScaleInPlace(float s, vector a) {
    a[0]*=s;
    a[1]*=s;
    a[2]*=s;
}

void vecCross(vector a, vector b, vector res) {
    res[0] = a[1] * b[2] - a[2] * b[1];
    res[1] = a[2] * b[0] - a[0] * b[2];
    res[2] = a[0] * b[1] - a[1] * b[0];
}

void vecCenter(vector a, vector b, vector res) {
    res[0] = (a[0] + b[0]) / 2.0;
    res[1] = (a[1] + b[1]) / 2.0;
    res[2] = (a[2] + b[2]) / 2.0;
}

typedef struct {
    int bondAtomIndices[2];
    int sat0AtomIndices[4];
    int nSat0;
    int sat3AtomIndices[4];
    int nSat3;
    double sumInertia;
    double sumLogInertia;
    float *wOmegaBuffer;
} t_rotBond;

typedef struct {
    /* constant data */
    int nAtoms;
    float resMass;
    int nRotBonds;
    float *atomMasses;
    float *atomMassesX3;
    t_rotBond *rotBonds;
    /* transient data */
    float *atomCrdMol;
    float *atomVelMol;
    float inertia[3];
    float prevAxes[9];
    float angMomLab[3];
    float angMomMol[3];
    float omegaMol[3];
    float omegaLab[3];
    float *COMposBuffer;
    float *COMposCurrent;
    float *COMposCorrRef;
    float *COMvelBuffer;
    float *COMvelCurrent;
    float *COMvelCorrRef;
    float *wOmegaBuffer;
    float *wOmegaCurrent;
    float *wOmegaCorrRef;
    /* accumulating data */
    double sumInertia[3];
    double sumLogInertia[3];
    double *totCorr;
    double *trCorr;
    double *rotCorr;
    double *rotBondCorr;
    int propCnt;
    int corrCnt;
    /* offset in MD info arrays */
    int offset;
} t_residue;

typedef struct {
    t_residue *residues;
    int nResidues;
    double *totCorr;
    double *trCorr;
    double *rotCorr;
    double *rotBondCorr;
    double inertia[3];
    double logInertia[3];
    double rotBondInertia;
    double logRotBondInertia;
    int postProcessed;
} t_residueList;

typedef struct {
    int frameSize;
    float *atomCrd;
    int bufferSize;
    float *atomVelBuffer;
    float *atomVelCurrent;
    float *atomVelCorrRef;
    int nCorr;
    int nAtomsSel;
} t_MDinfo;

typedef struct {
    int *indices;
    long long rank;
} t_dih;

int allocResidueList(t_residueList *residueList,int nRes,int nCorr) {
    int i;

    residueList->residues=(t_residue*)malloc(nRes*sizeof(t_residue));
    residueList->nResidues=nRes;
    residueList->totCorr=(double*)calloc(nCorr,sizeof(double));
    residueList->trCorr=(double*)calloc(nCorr*3,sizeof(double));
    residueList->rotCorr=(double*)calloc(nCorr*3,sizeof(double));
    residueList->rotBondCorr=(double*)calloc(nCorr,sizeof(double));
    for(i=0;i<3;i++) {
        residueList->inertia[i]=0.0;
        residueList->logInertia[i]=0.0;
    }
    residueList->rotBondInertia=0.0;
    residueList->logRotBondInertia=0.0;
    return 0;
}

int allocMDinfo(t_MDinfo *MDinfo, int nCorr, int nAtomsSel) {
    int i;

    MDinfo->frameSize=nAtomsSel*3;
    MDinfo->atomCrd=(float*)malloc(MDinfo->frameSize*sizeof(float));
    MDinfo->bufferSize=nCorr*MDinfo->frameSize;
    MDinfo->atomVelBuffer=(float*)malloc(MDinfo->bufferSize*sizeof(float));
    MDinfo->nCorr=nCorr;
    MDinfo->nAtomsSel=nAtomsSel;
    return 0;
}

int allocResidue(t_residue *res,int nAtoms,float *atomMasses,float resMass,int nAtomsSel,int nCorr) {
    int i,j;
    
    res->nAtoms=nAtoms;
    res->resMass=resMass;
    res->atomMasses=(float*)malloc(nAtoms*sizeof(float));
    res->atomMassesX3=(float*)malloc(nAtoms*3*sizeof(float));
    res->atomCrdMol=(float*)malloc(nAtoms*3*sizeof(float));
    res->atomVelMol=(float*)malloc(nAtoms*3*sizeof(float));
    for(i=0;i<nAtoms;i++) {
        res->atomMasses[i]=atomMasses[i];
        for(j=0;j<3;j++) {
            res->atomMassesX3[i*3+j]=atomMasses[i];
        }
    }
    for(i=0;i<3;i++) {
        res->sumInertia[i]=0.0;
        res->sumLogInertia[i]=0.0;
    }
    res->totCorr=(double*)calloc(nCorr, sizeof(double));
    res->COMposBuffer=(float*)calloc(nCorr*3, sizeof(float));
    res->COMvelBuffer=(float*)calloc(nCorr*3, sizeof(float));
    res->trCorr=(double*)calloc(nCorr*3, sizeof(double));
    res->wOmegaBuffer=(float*)calloc(nCorr*3, sizeof(float));
    res->rotCorr=(double*)calloc(nCorr*3, sizeof(double));
    for(i=0;i<9;i++) {
        res->prevAxes[i]=0.0;
    }
    res->prevAxes[0]=1.0;
    res->prevAxes[4]=1.0;
    res->prevAxes[8]=1.0;
    res->rotBondCorr=(double*)calloc(nCorr, sizeof(double));
    res->propCnt=0;
    res->corrCnt=0;
    return 0;
}

int setArrayIndexOffsets(t_residueList *residueList,t_MDinfo *MDinfo) {
    int i;
    int offset=0;

    for(i=0;i<residueList->nResidues;i++) {
        residueList->residues[i].offset=offset;
        offset+=residueList->residues[i].nAtoms*3;
    }
    return 0;
}

int computeTotalVACF(t_residue *res,t_MDinfo *MDinfo) {
    int i,j,k;
    float *m,*v1,*v2;

    m=res->atomMassesX3;
    v1=MDinfo->atomVelCorrRef+res->offset;
    for(i=0;i<MDinfo->nCorr;i++) {
        k = ((v1+(i*MDinfo->frameSize))-MDinfo->atomVelBuffer) % MDinfo->bufferSize;
        v2=MDinfo->atomVelBuffer+k;
        for(j=0;j<res->nAtoms*3;j++) {
            res->totCorr[i]+=m[j]*v1[j]*v2[j];
        }
    }
}

int computeCOM(int bufferShift,t_residue *res,t_MDinfo *MDinfo) {
    int i;
    float *m,*crd,*vel;

    m=res->atomMasses;
    res->COMposCurrent=res->COMposBuffer+bufferShift;
    res->COMvelCurrent=res->COMvelBuffer+bufferShift;
    crd=MDinfo->atomCrd+res->offset;
    vel=MDinfo->atomVelCurrent+res->offset;
    for(i=0;i<3;i++) {
        res->COMposCurrent[i]=0.0;
        res->COMvelCurrent[i]=0.0;
    }
    for(i=0;i<res->nAtoms;i++) {
        res->COMposCurrent[0]+=m[i]*crd[i*3+0];
        res->COMposCurrent[1]+=m[i]*crd[i*3+1];
        res->COMposCurrent[2]+=m[i]*crd[i*3+2];
        res->COMvelCurrent[0]+=m[i]*vel[i*3+0];
        res->COMvelCurrent[1]+=m[i]*vel[i*3+1];
        res->COMvelCurrent[2]+=m[i]*vel[i*3+2];
    }
    for(i=0;i<3;i++) {
        res->COMposCurrent[i]/=res->resMass;
        res->COMvelCurrent[i]/=res->resMass;
    }
    return 0;
}

int computeTransVACF(int bufferShift,t_residue *res,t_MDinfo *MDinfo) {
    int i,j,k;
    float mass;
    float *v1,*v2;

    mass=res->resMass;
    v1=res->COMvelBuffer+bufferShift;
    for(i=0;i<MDinfo->nCorr;i++) {
        k = ((v1+(i*3))-res->COMvelBuffer) % (MDinfo->nCorr*3);
        v2=res->COMvelBuffer+k;
        for(j=0;j<3;j++) {
            res->trCorr[i*3+j]+=mass*v1[j]*v2[j];
        }
    }
}

int subtractCOM(int bufferShift,t_residue *res,t_MDinfo *MDinfo) {
    int i;
    float *crdLab,*crd;
    float *velLab,*vel;
    float *COMpos,*COMvel;

    COMpos=res->COMposCurrent;
    COMvel=res->COMvelCurrent;
    crdLab=MDinfo->atomCrd+res->offset;
    velLab=MDinfo->atomVelCurrent+res->offset;
    crd=res->atomCrdMol;
    vel=res->atomVelMol;
    for(i=0;i<res->nAtoms;i++) {
        vecSub(crdLab,COMpos,crd);
        vecSub(velLab,COMvel,vel);
        crdLab+=3;
        velLab+=3;
        crd+=3;
        vel+=3;
    }
    return 0;
}

int computeRotation(int bufferShift,t_residue *res,t_MDinfo *MDinfo) {
    int i,j;
    float *m,*crd,*vel;
    float inertiaTensor[9];
    float axes[9];
    gsl_matrix *matrix;
    gsl_matrix *eVec;
    gsl_vector *eVal;
    gsl_eigen_symmv_workspace *ws;

    m=res->atomMasses;
    res->wOmegaCurrent=res->wOmegaBuffer+bufferShift;
    
    for(i=0;i<9;i++) {
        inertiaTensor[i]=0.0;
    }
    for(i=0;i<3;i++) {
        res->angMomLab[i]=0.0;
    }
    if(res->nAtoms>1) {
        crd=res->atomCrdMol;
        vel=res->atomVelMol;
        for(i=0;i<res->nAtoms;i++) {
            inertiaTensor[0]+=m[i]*(crd[1]*crd[1]+crd[2]*crd[2]);
            inertiaTensor[4]+=m[i]*(crd[0]*crd[0]+crd[2]*crd[2]);
            inertiaTensor[8]+=m[i]*(crd[0]*crd[0]+crd[1]*crd[1]);
            inertiaTensor[1]-=m[i]*crd[0]*crd[1];
            inertiaTensor[2]-=m[i]*crd[0]*crd[2];
            inertiaTensor[5]-=m[i]*crd[1]*crd[2];
            res->angMomLab[0]+=m[i]*(crd[1]*vel[2]-crd[2]*vel[1]);
            res->angMomLab[1]+=m[i]*(crd[2]*vel[0]-crd[0]*vel[2]);
            res->angMomLab[2]+=m[i]*(crd[0]*vel[1]-crd[1]*vel[0]);
            crd+=3;
            vel+=3;
        }
        inertiaTensor[3]=inertiaTensor[1];
        inertiaTensor[6]=inertiaTensor[2];
        inertiaTensor[7]=inertiaTensor[5];

        matrix=(gsl_matrix*)gsl_matrix_alloc(3,3);
        for(i=0;i<3;i++) {
            for(j=0;j<3;j++) {
                gsl_matrix_set(matrix,i,j,inertiaTensor[i*3+j]);
            }
        }
        eVal=(gsl_vector*)gsl_vector_alloc(3);
        eVec=(gsl_matrix*)gsl_matrix_alloc(3,3);
        /*allocating 'workspace'*/
        ws=(gsl_eigen_symmv_workspace*)gsl_eigen_symmv_alloc(12);
        /*compute eigenvalues and eigenvectors*/
        gsl_eigen_symmv(matrix,eVal,eVec,ws);
        /*sort eigenvalues and eigenvectors (largest eigenvalues first*/
        gsl_eigen_symmv_sort(eVal,eVec,GSL_EIGEN_SORT_VAL_DESC);

        for(i=0;i<3;i++) {
            res->inertia[i]=gsl_vector_get(eVal,i);
            res->sumInertia[i]+=res->inertia[i];
            if(res->inertia[i]>0.0) {
                res->sumLogInertia[i]+=log(res->inertia[i]);
            }
            for(j=0;j<3;j++) {
                axes[i*3+j]=gsl_matrix_get(eVec,j,i);
            }
            if(vecDot(&axes[i*3],&res->prevAxes[i*3])<0.0) {
                for(j=0;j<3;j++) {
                    axes[i*3+j]*=-1.0;
                }
            }
        }
        for(i=0;i<9;i++) {
            res->prevAxes[i]=axes[i];
        }
        for(i=0;i<3;i++) {
            res->angMomMol[i]=vecDot(&axes[i*3],res->angMomLab);
            if(res->inertia[i]>0.0) {
                res->omegaMol[i]=res->angMomMol[i]/res->inertia[i];
                res->wOmegaCurrent[i]=res->angMomMol[i]/sqrt(res->inertia[i]);
            } else {
                res->omegaMol[i]=0.0;
                res->wOmegaCurrent[i]=0.0;
            }
        }
        for(i=0;i<3;i++) {
            for(j=0;j<3;j++) {
                axes[i*3+j]=gsl_matrix_get(eVec,i,j);
            }
        }
        for(i=0;i<3;i++) {
            res->omegaLab[i]=vecDot(&axes[i*3],res->omegaMol);
        }
    }
    gsl_eigen_symmv_free(ws);
    gsl_vector_free(eVal);
    gsl_matrix_free(eVec);
    gsl_matrix_free(matrix);
    return 0;
}

int computeRotVACF(int bufferShift,t_residue *res,t_MDinfo *MDinfo) {
    int i,j,k;
    float *v1,*v2;

    v1=res->wOmegaBuffer+bufferShift;
    for(i=0;i<MDinfo->nCorr;i++) {
        k = ((v1+(i*3))-res->wOmegaBuffer) % (MDinfo->nCorr*3);
        v2=res->wOmegaBuffer+k;
        for(j=0;j<3;j++) {
            res->rotCorr[i*3+j]+=v1[j]*v2[j];
        }
    }
}

int subtractRot(t_residue *res) {
    int i;
    float *crd,*vel;
    float radVel[3];

    crd=res->atomCrdMol;
    vel=res->atomVelMol;
    for(i=0;i<res->nAtoms;i++) {
        vecCross(res->omegaLab, crd, radVel);
        vecSubInPlace(vel, radVel);
        crd+=3;
        vel+=3;
    }
    return 0;
}

int computeRotBonds(int bufferShift,t_residue *res) {
    int i,j;
    float b0[3], b1[3], sat[3], satVel[3], center[3], axis[3], proj[3], perp[3], mom[3], angMom[3];
    int nRotBonds=res->nRotBonds;
    t_rotBond *rotBond;
    float inertia1,inertia2;
    float reducedInertia;
    float angMom1,angMom2;
    float *crd,*vel,*masses;
    float *wOmegaCurrent;

    crd=res->atomCrdMol;
    vel=res->atomVelMol;
    masses=res->atomMasses;
     
    for(i=0;i<nRotBonds;i++) {
        inertia1=0.0;
        inertia2=0.0;
        angMom1=0.0;
        angMom2=0.0;
        rotBond=&res->rotBonds[i];
        wOmegaCurrent=rotBond->wOmegaBuffer+bufferShift;

        /*compute center of rotatable bond*/
        vecCopy(&crd[rotBond->bondAtomIndices[0]*3],b0);
        vecCopy(&crd[rotBond->bondAtomIndices[1]*3],b1);
        vecCenter(b0,b1,center);

        /*compute axis of rotation*/
        vecSub(b1,b0,axis);
        vecNormalizeInPlace(axis);

        /*compute inertia*/
        /*add inertia contributions from satellites with index 0 in dihedral*/
        for(j=0;j<rotBond->nSat0;j++) {
            vecCopy(&crd[rotBond->sat0AtomIndices[j]*3],sat);
	        vecCopy(&vel[rotBond->sat0AtomIndices[j]*3],satVel);
            /*coordinate of satellite relative to bond center*/
            vecSubInPlace(sat,center);
            /*project on rotational axis*/
            vecScale(vecDot(sat,axis),axis,proj);
            /*isolate perpendicular component*/
            vecSub(sat,proj,perp);
            /*moment of inertia*/
            inertia1+=masses[rotBond->sat0AtomIndices[j]]*vecDot(perp,perp);
            /*momentum*/
            vecScale(masses[rotBond->sat0AtomIndices[j]],satVel,mom);
            /*angular momentum with respect to axis*/
            vecCross(perp,mom,angMom);
            /*isolate angular momentum around axis*/
            angMom1+=vecDot(angMom,axis);
        }
        /*add inertia contributions from satellites with index 3 in dihedral*/
        for(j=0;j<rotBond->nSat3;j++) {
            vecCopy(&crd[rotBond->sat3AtomIndices[j]*3],sat);
	        vecCopy(&vel[rotBond->sat3AtomIndices[j]*3],satVel);
            /*coordinate of satellite relative to bond center*/
            vecSubInPlace(sat,center);
            /*project on rotational axis*/
            vecScale(vecDot(sat,axis),axis,proj);
            /*isolate perpendicular component*/
            vecSub(sat,proj,perp);
            /*moment of inertia*/
            inertia2+=masses[rotBond->sat3AtomIndices[j]]*vecDot(perp,perp);
            /*momentum*/
            vecScale(masses[rotBond->sat3AtomIndices[j]],satVel,mom);
            /*angular momentum with respect to axis*/
            vecCross(perp,mom,angMom);
            /*isolate angular momentum around axis*/
            angMom2+=vecDot(angMom,axis);
        }
        /*equivalent to reduced mass of vibrating bond in diatomic molecule*/
        reducedInertia=inertia1*inertia2/(inertia1+inertia2);
        /*accumulate reduced inertia for averaging*/
        rotBond->sumInertia+=reducedInertia;
        if(reducedInertia>0.0) {
            rotBond->sumLogInertia+=log(reducedInertia);
        }
	/*1) convert angular momenta into angular velocities*/
	/*2) compute difference (twist velocity of dihedral angle)*/
	/*3) multiply with square root of "reduced" inertia for weighted rotational velocity*/
        wOmegaCurrent[0]=sqrt(reducedInertia)*(angMom1/inertia1-angMom2/inertia2);
    }
    return 0;
}

int computeRotBondCorr(int bufferShift,t_residue *res,t_MDinfo *MDinfo) {
    int i,k,m;
    float *v1,*v2;

    for(m=0;m<res->nRotBonds;m++) {
        v1=res->rotBonds[m].wOmegaBuffer+bufferShift;
        for(i=0;i<MDinfo->nCorr;i++) {
            k = ((v1+i)-res->rotBonds[m].wOmegaBuffer) % MDinfo->nCorr;
            v2=res->rotBonds[m].wOmegaBuffer+k;
            res->rotBondCorr[i]+=v1[0]*v2[0];
        }
    }
    return 0;
}

int processStep(int tStep,t_MDinfo *MDinfo,float *crds,float *vels,t_residueList *residueList,int nCorr) {
    int i;
    int atomVelBufferShift;
    int resVecsBufferShiftCurrent;
    int resVecsBufferShiftCorrRef;
    int resNumBufferShiftCurrent;
    int resNumBufferShiftCorrRef;

    /* set MDinfo->atomVelCurrent for current time step */
    atomVelBufferShift=(tStep%MDinfo->nCorr)*MDinfo->nAtomsSel*3;
    MDinfo->atomVelCurrent=MDinfo->atomVelBuffer+atomVelBufferShift;
    /* set buffer shift for computation of residue property vectors*/
    resVecsBufferShiftCurrent=(tStep%MDinfo->nCorr)*3;
    resNumBufferShiftCurrent=(tStep%MDinfo->nCorr);

    #pragma omp parallel for
    for(i=0;i<3*MDinfo->nAtomsSel;i++) {
        MDinfo->atomCrd[i]=crds[i];
        MDinfo->atomVelCurrent[i]=vels[i];
    }

    if(tStep>=nCorr-1) {
        /* prep MDinfo->atomVelCorrRef for correlation functions*/
        atomVelBufferShift=((tStep+1)%MDinfo->nCorr)*MDinfo->nAtomsSel*3;
        MDinfo->atomVelCorrRef=MDinfo->atomVelBuffer+atomVelBufferShift;
        /* set buffer shift for correlation function of residue property vectors*/
        resVecsBufferShiftCorrRef=((tStep+1)%MDinfo->nCorr)*3;
        resNumBufferShiftCorrRef=((tStep+1)%MDinfo->nCorr);
    }
    
    #pragma omp parallel for
    for(i=0;i<residueList->nResidues;i++) {
        computeCOM(resVecsBufferShiftCurrent,&residueList->residues[i],MDinfo);
        subtractCOM(resVecsBufferShiftCurrent,&residueList->residues[i],MDinfo);
        computeRotation(resVecsBufferShiftCurrent,&residueList->residues[i],MDinfo);
        subtractRot(&residueList->residues[i]);
        computeRotBonds(resNumBufferShiftCurrent,&residueList->residues[i]);
        residueList->residues[i].propCnt++;
        if(tStep>=nCorr-1) {
            computeTotalVACF(&residueList->residues[i],MDinfo);
            computeTransVACF(resVecsBufferShiftCorrRef,&residueList->residues[i],MDinfo);
            computeRotVACF(resVecsBufferShiftCorrRef,&residueList->residues[i],MDinfo);
            computeRotBondCorr(resNumBufferShiftCorrRef,&residueList->residues[i],MDinfo);
            residueList->residues[i].corrCnt++;
        }
    }
    return 0;
}

int copyResidueListData(t_residueList *dst,t_residueList *src,int nCorr) {
    int i,j,k;

    for(j=0;j<nCorr;j++) {
        dst->totCorr[j]=src->totCorr[j];
        for(k=0;k<3;k++) {
            dst->trCorr[j*3+k]=src->trCorr[j*3+k];
            dst->rotCorr[j*3+k]=src->rotCorr[j*3+k];
        }
        dst->rotBondCorr[j]=src->rotBondCorr[j];
    }
    for(k=0;k<3;k++) {
        dst->inertia[k]=src->inertia[k];
        dst->logInertia[k]=src->logInertia[k];
    }
    dst->rotBondInertia=src->rotBondInertia;
    dst->logRotBondInertia=src->logRotBondInertia;

    for(i=0;i<src->nResidues;i++) {
        for(j=0;j<nCorr;j++) {
            dst->residues[i].totCorr[j]=src->residues[i].totCorr[j];
            for(k=0;k<3;k++) {
                dst->residues[i].trCorr[j*3+k]=src->residues[i].trCorr[j*3+k];
                dst->residues[i].rotCorr[j*3+k]=src->residues[i].rotCorr[j*3+k];
            }
            dst->residues[i].rotBondCorr[j]=src->residues[i].rotBondCorr[j];
        }
        for(k=0;k<3;k++) {
            dst->residues[i].sumInertia[k]=src->residues[i].sumInertia[k];
            dst->residues[i].sumLogInertia[k]=src->residues[i].sumLogInertia[k];
        }
        for(k=0;k<src->residues[i].nRotBonds;k++) {
            dst->residues[i].rotBonds[k].sumInertia=src->residues[i].rotBonds[k].sumInertia;
            dst->residues[i].rotBonds[k].sumLogInertia=src->residues[i].rotBonds[k].sumLogInertia;
        }
        dst->residues[i].corrCnt=src->residues[i].corrCnt;
        dst->residues[i].propCnt=src->residues[i].propCnt;
    }
    return 0;
}

int normResData(t_residue *res,int nCorr) {
    int j,k;
    int normFactor;

    normFactor=res->corrCnt;
    for(j=0;j<nCorr;j++) {
        res->totCorr[j]/=normFactor;
        for(k=0;k<3;k++) {
            res->trCorr[j*3+k]/=normFactor;
            res->rotCorr[j*3+k]/=normFactor;
        }
        res->rotBondCorr[j]/=normFactor;
    }
    normFactor=res->propCnt;
    for(k=0;k<3;k++) {
        res->sumInertia[k]/=normFactor;
        res->sumLogInertia[k]/=normFactor;
    }
    for(k=0;k<res->nRotBonds;k++) {
        res->rotBonds[k].sumInertia/=normFactor;
        res->rotBonds[k].sumLogInertia/=normFactor;
    }
    return 0;
}

int normResListData(t_residueList *residueList,int nCorr) {
    int j,k;
    int normFactor;

    normFactor=residueList->nResidues;
    for(j=0;j<nCorr;j++) {
        residueList->totCorr[j]/=normFactor;
        for(k=0;k<3;k++) {
            residueList->trCorr[j*3+k]/=normFactor;
            residueList->rotCorr[j*3+k]/=normFactor;
        }
        residueList->rotBondCorr[j]/=normFactor;
    }
    for(k=0;k<3;k++) {
        residueList->inertia[k]/=normFactor;
        residueList->logInertia[k]/=normFactor;
    }
    residueList->rotBondInertia/=normFactor;
    residueList->logRotBondInertia/=normFactor;
    return 0;
}

int sumData(t_residueList *residueList,int nCorr) {
    int i,j,k;
    int normFactor;

    for(i=0;i<residueList->nResidues;i++) {
        normResData(&residueList->residues[i],nCorr);
        for(j=0;j<nCorr;j++) {
            residueList->totCorr[j]+=residueList->residues[i].totCorr[j];
            for(k=0;k<3;k++) {
                residueList->trCorr[j*3+k]+=residueList->residues[i].trCorr[j*3+k];
                residueList->rotCorr[j*3+k]+=residueList->residues[i].rotCorr[j*3+k];
            }
            residueList->rotBondCorr[j]+=residueList->residues[i].rotBondCorr[j];
        }
        normFactor=residueList->residues[i].propCnt;
        for(k=0;k<3;k++) {
            residueList->inertia[k]+=residueList->residues[i].sumInertia[k];
            residueList->logInertia[k]+=residueList->residues[i].sumLogInertia[k];
        }
        for(k=0;k<residueList->residues[i].nRotBonds;k++) {
            residueList->rotBondInertia+=residueList->residues[i].rotBonds[k].sumInertia;
            residueList->logRotBondInertia+=residueList->residues[i].rotBonds[k].sumLogInertia;
        }
    }
    normResListData(residueList,nCorr);
    return 0;
}

/* qsort int comparison function */
int int_cmp(const void *a, const void *b)
{
        const int *ia = (const int *)a; /* casting pointer types  */
        const int *ib = (const int *)b;
        return *ia  - *ib;
        /* integer comparison: returns negative if b > a
        and positive if a > b */
}

/* qsort comparison function for dihedral index list */
int dih_cmp(const void *a, const void *b)
{
        const t_dih *ia = (const t_dih *)a; /* casting pointer types  */
        const t_dih *ib = (const t_dih *)b;

        long long diff = ia[0].rank  - ib[0].rank;

        if(diff<0) {
            return -1;
        } else if(diff>0) {
            return 1;
        } else {
            return 0;
        }
        /* integer comparison: returns negative if b > a
        and positive if a > b */
}

int getRotBonds(t_residueList *sets,int nSets,int *dihedAtomIndices,int nDih,int *resAtomIdxRange, int nCorr) {
    int i,j,k,l;
    int dih[4];
    int *tmp;
    int min,max;
    long long range;
    int cnt,idx;
    t_dih *dihList;
    t_rotBond *rotBonds;
    int flag;
    int nRotBonds;

    if(nDih>0) {
        dihList=(t_dih*)malloc(nDih*sizeof(t_dih));
        /*ensure order for dihedral list: [X a b Y] with a<b*/
        for(i=0;i<nDih;i++) {
            if(dihedAtomIndices[i*4+1] > dihedAtomIndices[i*4+2]) {
                for(k=0;k<4;k++) {
                    dih[k]=dihedAtomIndices[i*4+k];
                }
                for(k=0;k<4;k++) {
                    dihedAtomIndices[i*4+k]=dih[3-k];
                }
            }
        }
        
        /*find largest index in dihedral list*/
        min=dihedAtomIndices[0];
        max=-1;
        for(i=0;i<nDih*4;i++) {
            if(dihedAtomIndices[i]<min) {
                min=dihedAtomIndices[i];
            }
            if(dihedAtomIndices[i]>max) {
                max=dihedAtomIndices[i];
            }
        }
        range=(long long)(max-min+1);
        /*sort dihedral lists*/
        for(i=0;i<nDih;i++) {
            dihList[i].indices=&dihedAtomIndices[4*i];
            /*computing rank for sorting*/
            /*priority is on atom indices for rotatable bond*/
            dihList[i].rank =((long long)(dihList[i].indices[1]-min+1))*range*range*range;
            dihList[i].rank+=((long long)(dihList[i].indices[2]-min+1))*range*range;
            dihList[i].rank+=((long long)(dihList[i].indices[0]-min+1))*range;
            dihList[i].rank+=((long long)(dihList[i].indices[3]-min+1));
        }
        qsort(dihList,nDih,sizeof(t_dih),dih_cmp);
        tmp=(int*)malloc(nDih*sizeof(int));

        /*counting unique rotatable bonds: easy for sorted dihedrals*/
        i=0;
        j=dihList[i].indices[1];
        k=dihList[i].indices[2];
        cnt=1;
        for(i=1;i<nDih;i++) {
            if(dihList[i].indices[1]!=j || dihList[i].indices[2]!=k) {
                /*new unique rotBond detected*/
                j=dihList[i].indices[1];
                k=dihList[i].indices[2];
                cnt++;
            }
        }
        /*temporarily store all unique rotatable bonds in dihedral list in array 'rotBonds'*/
        rotBonds=(t_rotBond*)malloc(cnt*sizeof(t_rotBond));
        idx=0;
        i=0;
        j=dihList[i].indices[1];
        k=dihList[i].indices[2];
        rotBonds[idx].bondAtomIndices[0]=j;
        rotBonds[idx].bondAtomIndices[1]=k;
        rotBonds[idx].sat0AtomIndices[0]=dihList[i].indices[0];
        rotBonds[idx].nSat0=1;
        rotBonds[idx].sat3AtomIndices[0]=dihList[i].indices[3];
        rotBonds[idx].nSat3=1;
        for(i=1;i<nDih;i++) {
            if(dihList[i].indices[1]!=j || dihList[i].indices[2]!=k) {
                /*new unique rotBond detected*/
                j=dihList[i].indices[1];
                k=dihList[i].indices[2];
                idx++;
                rotBonds[idx].bondAtomIndices[0]=j;
                rotBonds[idx].bondAtomIndices[1]=k;
                rotBonds[idx].sat0AtomIndices[0]=dihList[i].indices[0];
                rotBonds[idx].nSat0=1;
                rotBonds[idx].sat3AtomIndices[0]=dihList[i].indices[3];
                rotBonds[idx].nSat3=1;
            } else {
                /*rotatable bond is the same*/
                /*test if the satellites are the same, otherwise add to list of satellites*/
                /*check satellite index 0*/
                flag=0;
                for(l=0;l<rotBonds[idx].nSat0;l++) {
                    if(dihList[i].indices[0]==rotBonds[idx].sat0AtomIndices[l]) {
                        flag=1;
                        break;
                    }
                }
                if(flag==0) {
                    if(rotBonds[idx].nSat0<4) {
                        rotBonds[idx].sat0AtomIndices[rotBonds[idx].nSat0]=dihList[i].indices[0];
                        rotBonds[idx].nSat0++;
                    } else {
                        printf("ERROR: found too many satellites for rotatable bond\n");
                        printf("       rotatable bond atom indices: %d %d\n",rotBonds[idx].bondAtomIndices[0],rotBonds[idx].bondAtomIndices[1]);
                        return 1;
                    }
                }
                /*check satellite index 3*/
                flag=0;
                for(l=0;l<rotBonds[idx].nSat3;l++) {
                    if(dihList[i].indices[3]==rotBonds[idx].sat3AtomIndices[l]) {
                        flag=1;
                        break;
                    }
                }
                if(flag==0) {
                    if(rotBonds[idx].nSat3<4) {
                        rotBonds[idx].sat3AtomIndices[rotBonds[idx].nSat3]=dihList[i].indices[3];
                        rotBonds[idx].nSat3++;
                    } else {
                        printf("ERROR: found too many satellites for rotatable bond\n");
                        printf("       rotatable bond atom indices: %d %d\n",rotBonds[idx].bondAtomIndices[0],rotBonds[idx].bondAtomIndices[1]);
                        return 1;
                    }
                }
            }
        }
        nRotBonds=idx+1;

        /*just for style, let's sort the satellite indices for each unique rotatable bond*/
        for(i=0;i<nRotBonds;i++) {
            qsort(rotBonds[i].sat0AtomIndices,rotBonds[i].nSat0,sizeof(int),int_cmp);
            qsort(rotBonds[i].sat3AtomIndices,rotBonds[i].nSat3,sizeof(int),int_cmp);
        }

        /*assign unique rotatable bonds to residues (one set of unique rotatable bonds per residue)*/
        /* loops over residues, i.e., sets*/
        // FILE *debug;
        // debug = fopen("/Users/mheyden/Documents/GitHub/HeydenLabASU/MDA-3D-2PT/debug.txt", "w");
        // fprintf(debug, "unique rotatable bonds\n");
        // for(i=0;i<nRotBonds;i++) {
        //     fprintf(debug, "rotBond %d\n", i+1);
        //     fprintf(debug, "bondAtomIndices: %d %d\n", rotBonds[i].bondAtomIndices[0], rotBonds[i].bondAtomIndices[1]);
        //     fprintf(debug, "sat0AtomIndices:");
        //     for(j=0;j<rotBonds[i].nSat0;j++) {
        //         fprintf(debug, " %d", rotBonds[i].sat0AtomIndices[j]);
        //     }
        //     fprintf(debug, "\n");
        //     fprintf(debug, "sat3AtomIndices:");
        //     for(j=0;j<rotBonds[i].nSat3;j++) {
        //         fprintf(debug, " %d", rotBonds[i].sat3AtomIndices[j]);
        //     }
        //     fprintf(debug, "\n");
        // }
        for(i=0;i<nSets;i++) {
            /*min max atom indices of residue i*/
            // fprintf(debug, "Residue %d\n", i+1);
            min=resAtomIdxRange[i*2+0];
            max=resAtomIdxRange[i*2+1];
            // fprintf(debug, "min: %d, max: %d\n", min, max);
            cnt=0;
            /*loop over unique rotatable bonds*/
            for(j=0;j<nRotBonds;j++) {
                flag=0;
                for(k=0;k<2;k++) {
                    l=rotBonds[j].bondAtomIndices[k];
                    if(l<min || l>max) {
                        flag=1;
                        break;
                    }
                }
                for(k=0;k<rotBonds[j].nSat0;k++) {
                    l=rotBonds[j].sat0AtomIndices[k];
                    if(l<min || l>max) {
                        flag=1;
                        break;
                    }
                }
                for(k=0;k<rotBonds[j].nSat3;k++) {
                    l=rotBonds[j].sat3AtomIndices[k];
                    if(l<min || l>max) {
                        flag=1;
                        break;
                    }
                }
                /*true if dihedral is part of residue i*/
                if(flag==0) {
                    tmp[cnt]=j;
                    // fprintf(debug, "rotBond %d %d\n",rotBonds[j].bondAtomIndices[0],rotBonds[j].bondAtomIndices[1]);
                    // fprintf(debug, "sat0AtomIndices:");
                    // for(k=0;k<rotBonds[j].nSat0;k++) {
                    //     fprintf(debug, " %d", rotBonds[j].sat0AtomIndices[k]);
                    // }
                    // fprintf(debug, "\n");
                    // fprintf(debug, "sat3AtomIndices:");
                    // for(k=0;k<rotBonds[j].nSat3;k++) {
                    //     fprintf(debug, " %d", rotBonds[j].sat3AtomIndices[k]);
                    // }
                    // fprintf(debug, "\n");
                    cnt++;
                }
            }
            /*allocate memory for rotatable bonds of residue i in sets->rotBondSets[i]*/
            sets->residues[i].nRotBonds=cnt;
            sets->residues[i].rotBonds=(t_rotBond*)malloc(cnt*sizeof(t_rotBond));
            /*copy rotatable bonds into for residue 'i' into sets->rotBondSets[i].rotBonds */
            for(j=0;j<cnt;j++) {
                l=tmp[j];
                for(k=0;k<2;k++) {
                    sets->residues[i].rotBonds[j].bondAtomIndices[k]=rotBonds[l].bondAtomIndices[k]-min;
                }
                for(k=0;k<rotBonds[l].nSat0;k++) {
                    sets->residues[i].rotBonds[j].sat0AtomIndices[k]=rotBonds[l].sat0AtomIndices[k]-min;
                }
                sets->residues[i].rotBonds[j].nSat0=rotBonds[tmp[j]].nSat0;
                for(k=0;k<rotBonds[l].nSat3;k++) {
                    sets->residues[i].rotBonds[j].sat3AtomIndices[k]=rotBonds[l].sat3AtomIndices[k]-min;
                }
                sets->residues[i].rotBonds[j].nSat3=rotBonds[tmp[j]].nSat3;
                sets->residues[i].rotBonds[j].sumInertia=0.0;
                sets->residues[i].rotBonds[j].sumLogInertia=0.0;
                sets->residues[i].rotBonds[j].wOmegaBuffer=(float*)malloc(nCorr*sizeof(float));
            }
        }
        // fclose(debug);
    } else {
        for(i=0;i<nSets;i++) {
            sets->residues[i].nRotBonds=0;
        }
    }
    return 0;
}
