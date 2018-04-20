#include "worker_code.h"

#ifdef SETTRAPFPE
#include <fenv.h>
#endif
#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <stddef.h>
#include <string.h>
#include <ctype.h>
#include <time.h>
#include <assert.h>
#include "mdl.h"
#include "master.h"
#include "outtype.h"
#include "smoothfcn.h"
#include "param.h"


/*DEBUG for FPE trapping... (not defined on most systems)*/
#ifdef SETTRAPFPE
#include <fpu_control.h>
#endif


MDL mdl;
MSR msr;
LCL *plcl;
PST pst;
PKD pkd;

FILE *fpLog = NULL;
FILE *fpLogTiming = NULL;
char achFile[256]; /*DEBUG use MAXPATHLEN here (& elsewhere)? -- DCR*/
double dTime;
double E=0,T=0,U=0,Eth=0,L[3]={0,0,0};
double dWMax=0,dIMax=0,dEMax=0,dMass=0,dMultiEff=0;
long lSec=0,lStart;
int i,j=0,iStep,bOutTime,iSec=0,iStop=0,nActive,nOutputList, OutputList[NUMOUTPUTS];
char achBaseMask[256];

// Replace file reading function, which also initializes particle store
double amuseInitStore(int nStore){// Replaces msrReadTipsy
//	struct dump h;
	struct inReadTipsy in;
    PST pst0 = msr->pst;
	LCL *plcl = pst0->plcl;
	double dTime,aTo,tTo,z;
	struct inSetParticleTypes intype;
	double sec,dsec;

	msr->N = 0;
	msr->nDark = 0;
	msr->nGas = 0;
	msr->nStar = 0;
	msr->nMaxOrderGas = msr->nGas - 1;
	msr->nMaxOrder = NIORDERGASBUFFER + msr->N - 1;  /* allow creating NIORDERGASBUFFER gas particles */
	msr->nMaxOrderDark = NIORDERGASBUFFER + msr->nGas + msr->nDark - 1;

	assert(msr->N == msr->nDark+msr->nGas+msr->nStar);

	dTime = 0.;
	in.dvFac = 1.0;

	msrSetStopStep(msr, dTime);

    // START replace msrOneNodeReadTipsy
    //nParts = malloc(msr->nThreads*sizeof(*nParts));

    // START replace pkdInitialize
    //PKD pkd;
    pkd = (PKD)malloc(sizeof(struct pkdContext));
    pkd->mdl = mdl;
    pkd->iOrder = msr->param.iOrder;
    pkd->idSelf = mdlSelf(mdl);
    pkd->nThreads = mdlThreads(mdl);
    pkd->nStore = nStore;
    pkd->nLocal = 0;
    pkd->nDark = msr->nDark;
    pkd->nGas = msr->nGas;
    pkd->nStar = msr->nStar;
    pkd->nMaxOrderGas = NIORDERGASBUFFER + msr->nGas - 1; 
    pkd->nMaxOrderDark = NIORDERGASBUFFER + msr->nGas + msr->nDark - 1;
    pkd->nRejects = 0;
    pkd->fPeriod[0] = msr->param.dxPeriod;
    pkd->fPeriod[1] = msr->param.dyPeriod;
    pkd->fPeriod[2] = msr->param.dzPeriod;
    pkd->dxInflow = msr->param.dxInflow;
    pkd->dxOutflow = msr->param.dxOutflow;
    pkd->sinkLog.nLog = 0;

    // Allocate memory for particles. Initially: no particles...
    pkd->pStore = mdlMalloc(pkd->mdl,(nStore+1)*sizeof(PARTICLE));
    pkd->kdNodes = NULL;
    pkd->piLeaf = NULL;
    pkd->kdTop = NULL;    

    /*
     ** Allocate initial interaction lists
     */
    pkd->nMaxPart = 500;
    pkd->nMaxCellSoft = 500;
    pkd->nMaxCellNewt = 500;
    pkd->nSqrtTmp = 500;
    pkd->ilp = malloc(pkd->nMaxPart*sizeof(ILP));
    mdlassert(pkd->mdl,pkd->ilp != NULL);
    pkd->ilcs = malloc(pkd->nMaxCellSoft*sizeof(ILCS));
    mdlassert(pkd->mdl,pkd->ilcs != NULL);
    pkd->ilcn = malloc(pkd->nMaxCellNewt*sizeof(ILCN));
    mdlassert(pkd->mdl,pkd->ilcn != NULL);
    pkd->sqrttmp = malloc(pkd->nSqrtTmp*sizeof(double));
    mdlassert(pkd->mdl,pkd->sqrttmp != NULL);
    pkd->d2a = malloc(pkd->nSqrtTmp*sizeof(double));
    mdlassert(pkd->mdl,pkd->d2a != NULL);
    /*
     ** Ewald stuff!
     */
    pkd->nMaxEwhLoop = 100;
    pkd->ewt = malloc(pkd->nMaxEwhLoop*sizeof(EWT));
    mdlassert(pkd->mdl,pkd->ewt != NULL);
    /*
     * Cooling
     */
#ifdef GASOLINE
#ifndef NOCOOLING
    pkd->Cool = CoolInit();
#endif  
#endif  

    plcl->pkd = pkd;

    // END replace pkdInitialize

    while(pst0->nLeaves > 1)
		pst0 = pst0->pstLower;

    // END replace msrOneNodeReadTipsy

    return dTime;
}

// Replacing param.c functions, since we want to avoid parsing parameter
// files and command line arguments.
void prmInitialize(PRM *pprm,void (*fcnLeader)(void),void (*fcnTrailer)(void))
{
	PRM prm;

	prm = (PRM)malloc(sizeof(struct prmContext));
	assert(prm != NULL);
	*pprm = prm;
	prm->pnHead = NULL;
	prm->fcnLeader = fcnLeader;
	prm->fcnTrailer = fcnTrailer;
	}

void prmFinish(PRM prm)
{
	PRM_NODE *pn,*pnKill;

	pn = prm->pnHead;
	while (pn) {
		pnKill = pn;
		pn = pn->pnNext;
		free(pnKill->pszName);
		if (pnKill->pszArg) free(pnKill->pszArg);
		if (pnKill->pszArgUsage) free(pnKill->pszArgUsage);
		free(pnKill);
		}
	free(prm);
	}

void prmAddParam(PRM prm,char *pszName,int iType,void *pValue,
				 int iSize,char *pszArg,char *pszArgUsage)
{
	PRM_NODE *pn,*pnTail;

	pn = (PRM_NODE *)malloc(sizeof(PRM_NODE));
	assert(pn != NULL);
	pn->pszName = (char *)malloc(strlen(pszName)+1);
	assert(pn->pszName != NULL);
	strcpy(pn->pszName,pszName);
	pn->iType = iType;
	pn->iSize = iSize;
	pn->bArg = 0;
	pn->bFile = 0;
	pn->pValue = pValue;
	if (pszArg) {
		pn->pszArg = (char *)malloc(strlen(pszArg)+1);
		assert(pn->pszArg != NULL);
		strcpy(pn->pszArg,pszArg);
		}
	else pn->pszArg = NULL;
	if (pszArgUsage) {
		pn->pszArgUsage = (char *)malloc(strlen(pszArgUsage)+1);
		assert(pn->pszArgUsage != NULL);
		strcpy(pn->pszArgUsage,pszArgUsage);
		}
	else pn->pszArgUsage = NULL;
	pn->pnNext = NULL;
	if (!prm->pnHead) prm->pnHead = pn;
	else {	
		pnTail = prm->pnHead;
		while (pnTail->pnNext) pnTail = pnTail->pnNext;
		pnTail->pnNext = pn;
		}
	}

void prmArgUsage(PRM prm)
{
	PRM_NODE *pn;

	if (prm->fcnLeader) (*prm->fcnLeader)();
	pn = prm->pnHead;
	while (pn) {
		if (pn->pszArg && pn->pszArgUsage) {
			if (pn->iType == 0) {
				printf("[+%s][-%s] %s\n",pn->pszArg,pn->pszArg,
					   pn->pszArgUsage);
				}
			else {
				printf("[-%s %s]\n",pn->pszArg,pn->pszArgUsage);
				}
			}
		pn = pn->pnNext;
		}
	if (prm->fcnTrailer) (*prm->fcnTrailer)();
	}

int prmParseParam(PRM prm,char *pszFile)
{
	PRM_NODE *pn;
	char achBuf[PRM_LINE_SIZE];
	char *p,*q,*pszCmd,t;
	int iLine,ret;
    
    // We must make sure msr->param.achInFile gets set to something
    // Doesn't matter what, we won't be reading it anyway
    // But let's be careful and point towards /dev/null
	p = strcpy(achBuf, "achInFile = /dev/null\n");
	iLine = 1;
	pszCmd = p;
    ++p;
	while (isalnum((int) *p)||strchr("_$",*p)) {
		++p;
		if (*p == 0) break;
		}
	t = *p;
	*p = 0;
	pn = prm->pnHead;
	while (pn) {
		if (!strcmp(pszCmd,pn->pszName)) break;
		pn = pn->pnNext;
		}
	*p = t;
	while (isspace((int) *p)) {
		++p;
		}
	++p;
	while (isspace((int) *p)) {
		++p;
		}
	ret = sscanf(p,"%[^\n#]",(char *)pn->pValue);
	pn->bFile = 1;
	return(1);
	}

int prmArgProc(PRM prm,int argc,char **argv){
	PRM_NODE *pn;
    char achBuf[PRM_LINE_SIZE];
    prmParseParam(prm,argv[argc-1]);

	return(1);
}

int prmArgSpecified(PRM prm,char *pszName)
{
	PRM_NODE *pn;
	
	pn = prm->pnHead;
	while (pn) {
		if (pn->pszArg)
			if (!strcmp(pn->pszName,pszName)) break;
		pn = pn->pnNext;
		}
	if (!pn) return(0);
	return(pn->bArg);
	}

int prmFileSpecified(PRM prm,char *pszName)
{
	PRM_NODE *pn;
	
	pn = prm->pnHead;	
	while (pn) {
		if (!strcmp(pn->pszName,pszName)) break;
		pn = pn->pnNext;
		}
	if (!pn) return(0);
	return(pn->bFile);
	}


int prmSpecified(PRM prm,char *pszName)
{
	return(prmArgSpecified(prm,pszName) || prmFileSpecified(prm,pszName));
	}

// from pkd.c
PARTICLE *p;

/* code to make gasoline core dump if there is a floating point exception 
 feenableexcept(FE_OVERFLOW | FE_DIVBYZERO | FE_INVALID);*/


#ifndef CCC
/* no stdout buffering */
//setbuf(stdout,(char *) NULL);
#endif

#ifdef __FAST_MATH__
assert(0); /* too dangerous! (gcc compile option) */
#endif

#ifdef SETTRAPFPE
feenableexcept(FE_INVALID|FE_DIVBYZERO|FE_OVERFLOW);
#endif

/*DEBUG following may work to explicitly trap FPEs -- also uncomment #include above */
/*
 #ifdef SETTRAPFPE
 {
 fpu_control_t cw = 0x1372;
 _FPU_SETCW(cw);
 }
 #endif
 */

int argc=0;
char *fake_argv[3];

void main_ch(MDL mdl)
{
    PST pst;
    LCL lcl;
    
    lcl.pszDataPath = (char *)getenv("PTOOLS_DATA_PATH");
    lcl.pkd = NULL;
    pstInitialize(&pst,mdl,&lcl);
    
    pstAddServices(pst,mdl);
    
    mdlHandler(mdl);
    
    pstFinish(pst);
}


int set_defaults(){
    return -1;
}


int initialize_code(){
    fake_argv[0] = "./gasoline_worker";
    fake_argv[1] = "/dev/null";
    fake_argv[2] = NULL;
    argc = 2;

    int nStore = 100; //FIXME: this is the maximum number of particles!
    
    lStart=time(0);
    //FIXME Remove argv things
    // mdlInitialize (MDL null/mpi variants) cares about two arguments:
    // "-sz", for nThreads,
    // "+d", for diagnostic output
    mdlInitialize(&mdl,fake_argv,main_ch);
    for(argc = 0; fake_argv[argc]; argc++); /* some MDLs can trash argv */
    msrInitialize(&msr,mdl,argc,fake_argv);
    
    dTime = amuseInitStore(nStore);

    printf("AMUSE: Gasoline initialized\n");

#ifdef GASOLINE
#ifndef NOCOOLING
    if (msr->param.iGasModel == GASMODEL_COOLING ||
        msr->param.bStarForm)
        msrInitCooling(msr);
        printf("AMUSE: Cooling initialized\n");
#ifdef OUTURBDRIVER
    msrInitouturb(msr, dTime);
    printf("AMUSE: outurb initialized\n");
#endif
    if(msr->param.bStarForm)
        msrInitStarLog(msr);
        printf("AMUSE: StarLog initialized\n");
#endif
    if(msr->param.bDoSinks && !msr->param.bBHSink)
        msrInitSinkLog(msr);
        printf("AMUSE: SinkLog initialized\n");
#endif
    if(msr->param.bRotatingBar)
        msrInitRotatingBar(msr, dTime);
        printf("AMUSE: RotatingBar initialized\n");
    msrInitStep(msr);
    printf("AMUSE: Step initialized\n");
#ifdef GLASS
    msrInitGlass(msr);
    printf("AMUSE: Glass initialized\n");
#endif
    dMass = msrMassCheck(msr,-1.0,"Initial");
    printf("AMUSE: Mass Checked\n");
    if (prmSpecified(msr->prm,"dSoft")) msrSetSoft(msr,msrSoft(msr));
    msrMassCheck(msr,dMass,"After msrSetSoft");
    printf("AMUSE: Mass Checked after SetSoft\n");
    
    msrSetSink(msr,dTime);
    printf("AMUSE: Sink Set\n");
#ifdef AGGS
    /* find and initialize any aggregates */
    msrAggsFind(msr);
    printf("AMUSE: Searched for Aggs\n");
    msrMassCheck(msr,dMass,"After msrAggsFind");
    printf("AMUSE: Mass Checked after AggsFind\n");
#endif
    /*
     ** If the simulation is periodic make sure to wrap all particles into
     ** the "unit" cell. Doing a drift of 0.0 will always take care of this.
     */
    msrDrift(msr,dTime,0.0); /* also finds initial overlaps for COLLISIONS */
    printf("AMUSE: Initial Drift\n");
    msrMassCheck(msr,dMass,"After initial msrDrift");
    printf("AMUSE: Mass Checked after initial Drift\n");
    
    /*
     ** Build tree, activating all particles first (just in case).
     */
    msrActiveType(msr,TYPE_ALL,TYPE_ACTIVE);
    msrDomainDecomp(msr,0,1);
    msrActiveType(msr,TYPE_ALL,TYPE_ACTIVE);
    msrInitAccel(msr);
    msrActiveType(msr,TYPE_ALL,TYPE_ACTIVE|TYPE_TREEACTIVE);
    msrUpdateSoft(msr,dTime);
    msrActiveType(msr,TYPE_ALL,TYPE_ACTIVE|TYPE_TREEACTIVE);
    msrBuildTree(msr,0,dMass,0);
    msrMassCheck(msr,dMass,"After msrBuildTree");
    if (msrDoGravity(msr)) {
        msrGravity(msr,0.0,msrDoSun(msr),&iSec,&dWMax,&dIMax,&dEMax,&nActive);
        msrMassCheck(msr,dMass,"After msrGravity");
        msrCalcEandL(msr,MSR_INIT_E,dTime,&E,&T,&U,&Eth,L);
        msrMassCheck(msr,dMass,"After msrCalcEandL");
        dMultiEff = 1.0;
    }
#ifdef GASOLINE
    msrInitSph(msr,dTime);
#endif
#ifdef INFLOWOUTFLOW
    if (msr->param.bInflowOutflow) msrModifyAccel(msr,dTime); /* zero acceleration of inflow/outflow particles */
#endif
    if (msr->param.bDoSinksAtStart) msrDoSinks(msr, dTime, 0.0, 0);
    
    return 0;
}

int commit_particles(){
    return 0;
}

int evolve_model(double time_end){        
    //LogTimingZeroCounters( msr );
    //FIXME Replacement for steps
    while(1) {
        iStep++;
        printf("dTime: %f  iStep: %i\n\n", dTime, iStep);
        if (dTime > time_end) {
            break;
        }
        if (msrKDK(msr)) {
            dMultiEff = 0.0;
            lSec = time(0);
            
            {
                msrTopStepKDK(msr,iStep-1,dTime,
                              msrDelta(msr),0,0,1,
                              &dMultiEff,&dWMax,&dIMax,
                              &dEMax,&iSec);
            }
            
            msrCoolVelocity(msr,dTime,dMass);	/* Supercooling if specified */
            msrMassCheck(msr,dMass,"After CoolVelocity in KDK");
            dTime += msrDelta(msr);
            
            lSec = time(0) - lSec;
        }
        else {
            lSec = time(0);
            msr->bDoneDomainDecomp = 0;
            msrTopStepDKD(msr,iStep-1,dTime,msrDelta(msr),&dMultiEff);
            msrRungStats(msr);
            msrCoolVelocity(msr,dTime,dMass); /* Supercooling if specified */
            msrMassCheck(msr,dMass,"After CoolVelocity in DKD");
            msrGrowMass(msr,dTime,msrDelta(msr)); /* Grow Masses if specified */
            dTime += msrDelta(msr);
            lSec = time(0) - lSec;
        }
    }
    return 0;
}

int cleanup_code(){
    msrFinish(msr);
    mdlFinish(mdl);
    return 0;
}

int get_mass(int index_of_the_particle, double * mass){
    return 0;
}

int get_time(double * time){
    *time = dTime;
    return 0;
}

int set_mass(int index_of_the_particle, double mass){
    return 0;
}

int get_index_of_first_particle(int * index_of_the_particle){
    return 0;
}

int get_total_radius(double * radius){
    return 0;
}

int get_potential_at_point(double eps, double x, double y, double z, 
                           double * phi, int npoints){
    return 0;
}

int get_total_mass(double * mass){
    return 0;
}

int set_eps2(double epsilon_squared){
    return 0;
}

int get_begin_time(double * time){
    //FIXME: no distinction between 'begin time' and 'time' yet...
    *time = dTime;
    return 0;
}

int get_eps2(double * epsilon_squared){
    return 0;
}

int get_index_of_next_particle(int index_of_the_particle, 
                               int * index_of_the_next_particle){
    return 0;
}

// int new_particle
//
// TODO: check msrReorder, ...
int new_dm_particle(int * index_of_the_particle, double mass, double x, 
                    double y, double z, double vx, double vy, double vz, 
                    double radius){
    msr->nDark += 1;
    msr->N += 1;
    msr->nMaxOrderDark = NIORDERGASBUFFER + msr->nGas + msr->nDark - 1;
    msr->nMaxOrder = NIORDERGASBUFFER + msr->N - 1;

    *index_of_the_particle = pkd->idSelf;
    
    PARTICLE p;
    p = pkd->pStore[pkd->idSelf];
    p.r[0] = x;
    p.r[1] = y;
    p.r[2] = z;
    p.v[0] = vx;
    p.v[1] = vy;
    p.v[2] = vz;
    p.fMass = mass;
    p.fSoft = radius;
    //p.fPot = phi;
    pkdNewParticle(pkd, p);
    return 0;
}

int new_sph_particle(int * index_of_the_particle, double mass, double x, 
                     double y, double z, double vx, double vy, double vz, 
                     double u, double h_smooth, double metals){
    msr->nGas += 1;
    msr->N += 1;
    msr->nMaxOrderGas = msr->nGas - 1;
    msr->nMaxOrderDark = NIORDERGASBUFFER + msr->nGas + msr->nDark - 1;
    msr->nMaxOrder = NIORDERGASBUFFER + msr->N - 1;

    *index_of_the_particle = msr->N;
    
    PARTICLE p;
    p = pkd->pStore[pkd->idSelf];
    TYPESet(&p, TYPE_GAS);
    p.r[0] = x;
    p.r[1] = y;
    p.r[2] = z;
    p.v[0] = vx;
    p.v[1] = vy;
    p.v[2] = vz;
    p.fMass = mass;
    p.fSoft = h_smooth;
    //p.fPot = phi;
    p.u = u;
    p.uPred = u;
    p.fMetals = metals;
    //p.fMetalsPred = metals;
    pkdNewParticle(pkd, p);
    return 0;
}

int new_star_particle(int * index_of_the_particle, double mass, double x, 
                      double y, double z, double vx, double vy, double vz,
                      double tform, double radius, double metals){

    msr->nStar += 1;
    msr->N += 1;
    msr->nMaxOrder = NIORDERGASBUFFER + msr->N - 1;
    
    *index_of_the_particle = msr->N;
    
    PARTICLE p;
    p = pkd->pStore[pkd->idSelf];
    TYPESet(&p, TYPE_STAR); // should set to TYPE_SINK for black holes
    p.r[0] = x;
    p.r[1] = y;
    p.r[2] = z;
    p.v[0] = vx;
    p.v[1] = vy;
    p.v[2] = vz;
    p.fMass = mass;
    p.fMassForm = mass;
    p.fSoft = radius;
    p.fSoft0 = radius;
    p.fTimeForm = tform;
    p.fBallMax = 0.0;
    p.fMetals = metals;
    p.fNSNtot = 0.0;

    pkdNewParticle(pkd, p);
    return 0;
}

int delete_particle(int index_of_the_particle){
    // Check if particle is dm, star or sph
    // make sure index_of_the_particle is adjusted for all following particles
    //
    // pkdDeleteParticle(pkd,&pkd->pStore[pDel->iIndex]);
    msr->N -= 1;
    return -1;
}

int get_potential(int index_of_the_particle, double * potential){
    return 0;
}

int synchronize_model(){
    return 0;
}

int set_state(int index_of_the_particle, double mass, double x, double y, 
              double z, double vx, double vy, double vz, double radius){
    return 0;
}

int get_state(int index_of_the_particle, double * mass, double * x, 
              double * y, double * z, double * vx, double * vy, double * vz,
              double * radius){
    int id = msr->pMap[index_of_the_particle];
    p = &pkd->pStore[id];
    * mass = p->fMass;
    * x = p->r[0];
    * y = p->r[1];
    * z = p->r[2];
    * vx = p->v[0];
    * vy = p->v[1];
    * vz = p->v[2];
    * radius = p->fSoft;
    return 0;
}

int set_state_sph(int index_of_the_particle, double mass, double x, double y, 
                  double z, double vx, double vy, double vz, double u, 
                  double h_smooth, double metals){
    return 0;
}

int get_state_sph(int index_of_the_particle, double * mass, double * x, 
                  double * y, double * z, double * vx, double * vy, double * vz,
                  double * u, double * h_smooth, double * metals){
    int id = msr->pMap[index_of_the_particle];
    p = &pkd->pStore[id];
    * mass = p->fMass;
    * x = p->r[0];
    * y = p->r[1];
    * z = p->r[2];
    * vx = p->v[0];
    * vy = p->v[1];
    * vz = p->v[2];
    * u = p->u;
    * h_smooth = p->fSoft;
    * metals = p->fMetals;
    return 0;
}

int set_state_star(int index_of_the_particle, double mass, double x, double y, 
                   double z, double vx, double vy, double vz, double tform,
                   double radius, double metals){
    return 0;
}

int get_state_star(int index_of_the_particle, double * mass, double * x, 
                   double * y, double * z, double * vx, double * vy, double * vz,
                   double * tform, double * radius, double * metals){
    int id = msr->pMap[index_of_the_particle];
    p = &pkd->pStore[id];
    * mass = p->fMass;
    * x = p->r[0];
    * y = p->r[1];
    * z = p->r[2];
    * vx = p->v[0];
    * vy = p->v[1];
    * vz = p->v[2];
    * tform = p->fTimeForm;
    * radius = p->fSoft;
    * metals = p->fMetals;
    return 0;
}

int get_time_step(double * time_step){
    *time_step = msr->param.dDelta;
    return 0;
}

int set_time_step(double time_step){
    msr->param.dDelta = time_step;
    return 0;
}

int recommit_particles(){
    return 0;
}

int get_kinetic_energy(double * kinetic_energy){
    return 0;
}

int get_number_of_particles(int * number_of_particles){
    *number_of_particles = msr->N;
    return 0;
}

int set_acceleration(int index_of_the_particle, double ax, double ay, 
                     double az){
    return 0;
}

int get_center_of_mass_position(double * x, double * y, double * z){
    return 0;
}

int get_center_of_mass_velocity(double * vx, double * vy, double * vz){
    return 0;
}

int get_radius(int index_of_the_particle, double * radius){
    return 0;
}

int get_metallicity(int index_of_the_particle, double * metals){
    return 0;
}

int get_star_tform(int index_of_the_particle, double * tform){
    return 0;
}

int set_begin_time(double time){
    return 0;
}

int set_radius(int index_of_the_particle, double radius){
    return 0;
}

int set_metallicity(int index_of_the_particle, double metals){
    return 0;
}

int set_star_tform(int index_of_the_particle, double tform){
    return 0;
}

int recommit_parameters(){
    return 0;
}

int get_potential_energy(double * potential_energy){
    return 0;
}

int get_gravity_at_point(double eps, double x, double y, double z, 
                         double * ax, double * ay, double * az, int npoints){
    return 0;
}

int get_velocity(int index_of_the_particle, double * vx, double * vy, double * vz){
    return 0;
}

int get_internal_energy(int index_of_the_particle, double * u){
    return 0;
}

int get_position(int index_of_the_particle, double * x, double * y, double * z){
    PST pst0;
    LCL *plcl;
    int i = 1;
    pst0 = msr->pst;
    plcl = pst0->plcl;
    *x = plcl->pkd->pStore[i].r[0];
    *y = plcl->pkd->pStore[i].r[1];
    *z = plcl->pkd->pStore[i].r[2];
    return 0;
}

int set_position(int index_of_the_particle, double x, double y, double z){
    return 0;
}

int get_acceleration(int index_of_the_particle, double * ax, double * ay, 
                     double * az){
    return 0;
}

int commit_parameters(){
    return 0;
}

int set_velocity(int index_of_the_particle, double vx, double vy, double vz){
    return 0;
}

int set_internal_energy(int index_of_the_particle, double u){
    return 0;
}

int get_alpha(double * alpha){
    *alpha = msr->param.dConstAlpha;
    return 0;
}

int set_alpha(double alpha){
    msr->param.dConstAlpha = alpha;
    return 0;
}

int get_beta(double * beta){
    *beta = msr->param.dConstBeta;
    return 0;
}

int set_beta(double beta){
    msr->param.dConstBeta = beta;
    return 0;
}

int get_eta(double * eta){
    *eta = msr->param.dEta;
    return 0;
}

int set_eta(double eta){
    msr->param.dEta = eta;
    return 0;
}

int get_eta_courant(double * eta_courant){
    *eta_courant = msr->param.dEtaCourant;
    return 0;
}

int set_eta_courant(double eta_courant){
    msr->param.dEtaCourant = eta_courant;
    return 0;
}

int get_theta(double * theta){
    *theta = msr->param.dTheta;
    return 0;
}

int set_theta(double theta){
    msr->param.dTheta = theta;
    return 0;
}

int get_theta_2(double * theta_2){
    *theta_2 = msr->param.dTheta2;
    return 0;
}

int set_theta_2(double theta_2){
    msr->param.dTheta2 = theta_2;
    return 0;
}

int get_nsmooth(int * nsmooth){
    *nsmooth = msr->param.nSmooth;
    return 0;
}

int set_nsmooth(int nsmooth){
    msr->param.nSmooth = nsmooth;
    return 0;
}
