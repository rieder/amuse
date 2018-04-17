#include "worker_code.h"

#ifdef SETTRAPFPE
#include <fenv.h>
#endif
#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <assert.h>
#include "mdl.h"
#include "master.h"
#include "outtype.h"
#include "smoothfcn.h"

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
    fake_argv[1] = "param.dat";
    fake_argv[2] = NULL;
    argc = 2;
    
    lStart=time(0);
    //FIXME Remove argv things
    mdlInitialize(&mdl,fake_argv,main_ch);
    for(argc = 0; fake_argv[argc]; argc++); /* some MDLs can trash argv */
    //FIXME: Replace this function (in master.c) by individual setters
    msrInitialize(&msr,mdl,argc,fake_argv);
    
    /*
     ** Read in the binary file, this may set the number of timesteps or
     ** the size of the timestep when the zto parameter is used.
     */
    dTime = msrReadTipsy(msr);
    
#ifdef GASOLINE
#ifndef NOCOOLING
    if (msr->param.iGasModel == GASMODEL_COOLING ||
        msr->param.bStarForm)
        msrInitCooling(msr);
#ifdef OUTURBDRIVER
    msrInitouturb(msr, dTime);
#endif
    if(msr->param.bStarForm)
        msrInitStarLog(msr);
#endif
    if(msr->param.bDoSinks && !msr->param.bBHSink)
        msrInitSinkLog(msr);
#endif
    if(msr->param.bRotatingBar)
        msrInitRotatingBar(msr, dTime);
    msrInitStep(msr);
#ifdef GLASS
    msrInitGlass(msr);
#endif
    dMass = msrMassCheck(msr,-1.0,"Initial");
    if (prmSpecified(msr->prm,"dSoft")) msrSetSoft(msr,msrSoft(msr));
    msrMassCheck(msr,dMass,"After msrSetSoft");
    
    msrSetSink(msr,dTime);
#ifdef AGGS
    /* find and initialize any aggregates */
    msrAggsFind(msr);
    msrMassCheck(msr,dMass,"After msrAggsFind");
#endif
    /*
     ** If the simulation is periodic make sure to wrap all particles into
     ** the "unit" cell. Doing a drift of 0.0 will always take care of this.
     */
    msrDrift(msr,dTime,0.0); /* also finds initial overlaps for COLLISIONS */
    msrMassCheck(msr,dMass,"After initial msrDrift");
    
    //    if (msrSteps(msr) > 0) {
    if (msrComove(msr)) {
        msrSwitchTheta(msr,dTime);
    }
    /*
     ** Now we have all the parameters for the simulation we can make a
     ** log file entry.
     */
    //if (msrLogInterval(msr)) {
    //    //FIXME Remove logging to file
    //    //sprintf(achFile,"%s.log",msrOutName(msr));
    //    fpLog = fopen(achFile,"w");
    //    assert(fpLog != NULL);
    //    setbuf(fpLog,(char *) NULL); /* no buffering */
    //    /*
    //     ** Include a comment at the start of the log file showing the
    //     ** command line options.
    //     */
    //    //FIXME Change logging from file to stderr
    //    fprintf(fpLog,"# ");
    //    //FIXME Remove argv things
    //    for (i=0;i<argc;++i) fprintf(fpLog,"%s ",fake_argv[i]);
    //    fprintf(fpLog,"\n");
    //    msrLogHeader(msr,fpLog);
    //    /* Timing data, if requested */
    //    fpLogTiming = LogTimingInit( msr, "w" );
    //}
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
        //if (msrLogInterval(msr)) {
        //    //FIXME Change logging from file to stderr
        //    (void) fprintf(fpLog,"%.4e %.4e %.6e %.4e %.4e %.4e %.6e %.6e %.6e "
        //                   "%i %.4e %.4e %.4e %.4e\n",dTime,
        //                   1.0/csmTime2Exp(msr->param.csm,dTime)-1.0,
        //                   E,T,U,Eth,L[0],L[1],L[2],iSec,dWMax,dIMax,dEMax,
        //                   dMultiEff);
        //}
        /* LogTimingOutput( msr, fpLogTiming, dTime, 0 ); */
    }
#ifdef GASOLINE
    msrInitSph(msr,dTime);
#endif
#ifdef INFLOWOUTFLOW
    if (msr->param.bInflowOutflow) msrModifyAccel(msr,dTime); /* zero acceleration of inflow/outflow particles */
#endif
    if (msr->param.bDoSinksAtStart) msrDoSinks(msr, dTime, 0.0, 0);
    /*
     ** Dump Frame Initialization
     */
    //msrDumpFrameInit( msr, dTime, 1.0*msr->param.iStartStep, 0);
    
    return 0;
}

int evolve_model(double time_end){        
    //LogTimingZeroCounters( msr );
    //FIXME Replacement for steps
    while(1) {
        iStep++;
        printf("dTime: %f  iStep: %i\n\n", dTime, iStep);
        //dTime = (double)iStep;
        //printf("recalculated dTime: %f  iStep: %i\n\n", dTime, iStep);
        if (dTime > time_end) {
            break;
        }
        //for (iStep=msr->param.iStartStep+1;iStep<=msr->param.iStopStep;++iStep) {
        if (msrComove(msr)) {
            msrSwitchTheta(msr,dTime);
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
            
            /* msrRungStats(msr); This is useless */
            msrCoolVelocity(msr,dTime,dMass);	/* Supercooling if specified */
            msrMassCheck(msr,dMass,"After CoolVelocity in KDK");
            dTime += msrDelta(msr);
            //if(iStep%msr->param.iOrbitOutInterval == 0) {
            //    msrOutputBlackHoles(msr, dTime);
            //}
            
            lSec = time(0) - lSec;
            /*
             ** Output a log file line if requested.
             ** Note: no extra gravity calculation required.
             */
            //if (msrLogInterval(msr) && iStep%msrLogInterval(msr) == 0) {
            //    msrCalcEandL(msr,MSR_STEP_E,dTime,&E,&T,&U,&Eth,L);
            //    msrMassCheck(msr,dMass,"After msrCalcEandL in KDK");
            //    //FIXME Change logging from file to stderr
            //    (void) fprintf(fpLog,"%.4e %.4e %.6e %.4e %.4e %.4e %.6e %.6e %.6e "
            //                   "%li %.4e %.4e %.4e %.4e\n",dTime,
            //                   1.0/csmTime2Exp(msr->param.csm,dTime)-1.0,
            //                   E,T,U,Eth,L[0],L[1],L[2],lSec,dWMax,dIMax,dEMax,
            //                   dMultiEff);
            //}
            //LogTimingOutput( msr, fpLogTiming, dTime, 0 );
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
            //if (msrLogInterval(msr) && iStep%msrLogInterval(msr) == 0) {
            //    /*
            //     ** Output a log file line.
            //     ** Reactivate all particles.
            //     */
            //    if (msrDoGravity(msr)) {
            //        msrActiveType(msr,TYPE_ALL,TYPE_ACTIVE|TYPE_TREEACTIVE);
            //        msrDomainDecomp(msr,0,1);
            //        msrActiveType(msr,TYPE_ALL,TYPE_ACTIVE|TYPE_TREEACTIVE);
            //        msrUpdateSoft(msr,dTime);
            //        msrActiveType(msr,TYPE_ALL,TYPE_ACTIVE|TYPE_TREEACTIVE);
            //        msrBuildTree(msr,0,dMass,0);
            //        msrMassCheck(msr,dMass,"After msrBuildTree in DKD-log");
            //        msrInitAccel(msr);
            //        msrGravity(msr,iStep,msrDoSun(msr),&iSec,&dWMax,&dIMax,&dEMax,&nActive);
            //        msrMassCheck(msr,dMass,"After msrGravity in DKD-log");
            //    }
            //    msrCalcEandL(msr,MSR_STEP_E,dTime,&E,&T,&U,&Eth,L);
            //    msrMassCheck(msr,dMass,"After msrCalcEandL in DKD-log");
            //    //FIXME Change logging from file to stderr
            //    (void) fprintf(fpLog,"%.4e %.4e %.6e %.4e %.4e %.4e %.6e %.6e %.6e "
            //                   "%li %.4e %.4e %.4e %.4e\n",dTime,
            //                   1.0/csmTime2Exp(msr->param.csm,dTime)-1.0,
            //                   E,T,U,Eth,L[0],L[1],L[2],time(0)-lSec,dWMax,dIMax,dEMax,
            //                   dMultiEff);
            //    
            //}
            //LogTimingOutput( msr, fpLogTiming, dTime, 0 );
            lSec = time(0) - lSec;
        }
        /*
         ** Check for user interrupt.
         */
        //iStop = msrCheckForStop(msr);
        /*
         ** Output if 1) we've hit an output time
         **           2) We are stopping
         **           3) we're at an output interval
         */
        //if (msr->param.iTreeZipStep && (iStep % msr->param.iTreeZipStep)==0) msrTreeZip(msr,iStep);
        
        //if ((bOutTime=msrOutTime(msr,dTime)) || iStep == msr->param.iStopStep || iStop ||
        //    (msrOutInterval(msr) > 0 && iStep%msrOutInterval(msr) == 0)  ||
        //    (msr->param.iOutMinorInterval && (iStep%msr->param.iOutMinorInterval == 0))) {
        //    int bDensitySmooth;
        //    if (msr->nGas && !msr->param.bKDK) {
        //        msrActiveType(msr,TYPE_GAS,TYPE_ACTIVE|TYPE_TREEACTIVE|TYPE_SMOOTHACTIVE);
        //        msrBuildTree(msr,1,-1.0,1);
        //        msrSmooth(msr,dTime,SMX_DENSITY,1);
        //    }
        //    bDensitySmooth = msrDoDensity(msr) || msr->param.bDohOutput;
        //    msrSelectOutputList(msr, &nOutputList, OutputList, iStep, bOutTime, &bDensitySmooth);
        //    
        //    if (bDensitySmooth) {
        //        msrActiveType(msr,TYPE_ALL,TYPE_ACTIVE|TYPE_TREEACTIVE|TYPE_SMOOTHACTIVE);
        //        msrDomainDecomp(msr,0,1);
        //        msrBuildTree(msr,0,dMass,1);
        //        msrSmooth(msr,dTime,SMX_DENSITYTMP,1);
        //        if (!msr->param.bNoReOrder) msrReorder(msr);
        //    }
        //    if (!msr->param.bNoReOrder) {
        //        msrReorder(msr);
        //        msrMassCheck(msr,dMass,"After msrReorder in OutTime");
        //    }
        //    //FIXME Remove logging to file
        //    //sprintf(achFile,msr->param.achDigitMask,msrOutName(msr),iStep);
        //    //msrWriteOutputs(msr, achFile, OutputList, nOutputList, dTime);
        //    //msrFlushStarLog(msr);
        //    //msrFlushSinkLog(msr);
        //    /*
        //     ** Don't allow duplicate outputs.
        //     */
        //    while (msrOutTime(msr,dTime));
        //}
        //if (!iStop && msr->param.iWallRunTime > 0) {
        //    if (msr->param.iWallRunTime*60 - (time(0)-lStart) < ((int) (lSec*1.5)) ) {
        //        printf("RunTime limit exceeded.  Writing checkpoint and exiting.\n");
        //        printf("    iWallRunTime(sec): %d   Time running: %ld   Last step: %ld\n",
        //               msr->param.iWallRunTime*60,time(0)-lStart,lSec);
        //        iStop = 1;
        //    }
        //}
        //if (iStop || iStep == msr->param.iStopStep ||
        //    (msrCheckInterval(msr) && iStep%msrCheckInterval(msr) == 0)) {
        //    /*
        //     ** Write a checkpoint.
        //     */
        //    msrFlushStarLog(msr);
        //    msrFlushSinkLog(msr);
        //    msrWriteCheck(msr,dTime,iStep);
        //    msrMassCheck(msr,dMass,"After msrWriteCheck");
        //Restart:
        //    ;
        //}
        //if (iStop) break;
    }
    //if (msrLogInterval(msr)) {
    //    (void) fclose(fpLog);
    //    LogTimingFinish( msr, fpLogTiming, dTime );
    //}
    //if (msr->param.bVStart) printf("Integration complete\n");
    //    }
    //    else {
    //        /* Do DiagnosticOutput */
    //        struct inInitDt in;
    //        msrActiveType(msr,TYPE_ALL,TYPE_ACTIVE|TYPE_TREEACTIVE|TYPE_SMOOTHACTIVE);
    //        
    //        in.dDelta = 1e37; /* large number */
    //        pstInitDt(msr->pst,&in,sizeof(in),NULL,NULL);
    //        msrInitAccel(msr);
    //        //puts("Initialized Accel and dt\n");
    //        //FIXME Remove logging to file
    //        //sprintf(achFile,"%s",msrOutName(msr));
    //        
    //        if (msrRestart(msr)) {
    //            msrReorder(msr);
    //            //FIXME Remove logging to file
    //            //sprintf(achFile,"%s",msrOutName(msr));
    //            //FIXME Writing output here
    //            //msrWriteTipsy(msr,achFile,dTime);
    //        }
    //        
    //        if (msrDoGravity(msr)) {
    //            msrActiveType(msr,TYPE_ALL,TYPE_ACTIVE|TYPE_TREEACTIVE|TYPE_SMOOTHACTIVE );
    //            msrDomainDecomp(msr,0,1);
    //            msrUpdateSoft(msr,dTime);
    //            msrBuildTree(msr,0,dMass,0);
    //            msrMassCheck(msr,dMass,"After msrBuildTree in OutSingle Gravity");
    //            msrGravity(msr,0.0,msrDoSun(msr),&iSec,&dWMax,&dIMax,&dEMax,&nActive);
    //            msrMassCheck(msr,dMass,"After msrGravity in OutSingle Gravity");
    //            msrReorder(msr);
    //            msrMassCheck(msr,dMass,"After msrReorder in OutSingle Gravity");
    //            nOutputList = 0;
    //            OutputList[nOutputList++]=OUT_ACCELG_VECTOR;
    //            OutputList[nOutputList++]=OUT_POT_ARRAY;
    //            msrWriteOutputs(msr, achFile, OutputList, nOutputList, dTime);
    //            msrMassCheck(msr,dMass,"After msrOutArray in OutSingle Gravity");
    //        }
    //#ifdef GASOLINE
    //        if (msr->nGas > 0) {
    //            msrInitSph(msr,dTime);
    //            msrCreateGasStepZeroOutputList(msr, &nOutputList,OutputList);
    //            OutputList[(nOutputList)++]=OUT_METALSDOT_ARRAY;
    //#ifdef STARFORM
    //            OutputList[(nOutputList)++]=OUT_OXYGENMASSFRACDOT_ARRAY;
    //            OutputList[(nOutputList)++]=OUT_IRONMASSFRACDOT_ARRAY;
    //#endif
    //#ifndef NOCOOLING
    //            
    //            if (msr->param.bGasCooling) {
    //                OutputList[nOutputList++]=OUT_COOL_EDOT_ARRAY;
    //                OutputList[nOutputList++]=OUT_COOL_COOLING_ARRAY;
    //                OutputList[nOutputList++]=OUT_COOL_HEATING_ARRAY;
    //            }
    //#endif
    //            if (msr->param.bSphStep) {
    //                //FIXME Change logging from file to stderr
    //                //fprintf(stdout,"Adding SphStep dt\n");
    //                msrSphStep(msr,dTime,0);
    //            }
    //            msrReorder(msr);
    //            //msrWriteOutputs(msr, achFile, OutputList, nOutputList, dTime);
    //            //msrFlushStarLog(msr);
    //            //msrFlushSinkLog(msr);
    //        }
    //#endif /* GASOLINE */
    //        /*
    //         ** Build tree, activating all particles first (just in case).
    //         */
    //        if (msrDoDensity(msr) || msr->param.bDensityStep) {
    //            msrActiveType(msr,TYPE_ALL,TYPE_ACTIVE|TYPE_TREEACTIVE|TYPE_SMOOTHACTIVE);
    //            msrDomainDecomp(msr,0,1);
    //            msrBuildTree(msr,0,-1.0,1);
    //            msrMassCheck(msr,dMass,"After msrBuildTree in OutSingle Density");
    //            //printf("Calculating DENSITYTMP\n");
    //            /* Note: this is a gather-scatter density -- different to SPH DENDVDX
    //             (same if we used msrSmooth(msr,dTime,SMX_DENSITYTMP,0) */
    //            msrSmooth(msr,dTime,SMX_DENSITYTMP,1);
    //            msrMassCheck(msr,dMass,"After msrSmooth in OutSingle Density");
    //        }
    //        if (msrDoDensity(msr)) {
    //            msrReorder(msr);
    //            msrMassCheck(msr,dMass,"After msrReorder in OutSingle Density");
    //            //nOutputList = 0;
    //            //OutputList[nOutputList++]=OUT_DENSITY_ARRAY;
    //            //msrWriteOutputs(msr, achFile, OutputList, nOutputList, dTime);
    //            /*	sprintf(achFile,"%s.den",msrOutName(msr));
    //             msrReorder(msr);
    //             msrOutArray(msr,achFile,OUT_DENSITY_ARRAY);*/
    //            msrMassCheck(msr,dMass,"After msrOutArray in OutSingle Density");
    //            
    //        }
    //        
    //        if (msrDoGravity(msr)) {
    //            if (msr->param.bGravStep) {
    //                //FIXME Change logging from file to stderr
    //                //fprintf(stdout,"Adding GravStep dt\n");
    //                msrGravStep(msr,dTime);
    //            }
    //            if (msr->param.bAccelStep) {
    //                //FIXME Change logging from file to stderr
    //                //fprintf(stdout,"Adding AccelStep dt\n");
    //                msrAccelStep(msr,dTime);
    //            }
    //        }
    //        
    //        if (msr->param.bDensityStep) {
    //            //FIXME Change logging from file to stderr
    //            //fprintf(stdout,"Adding DensStep dt\n");
    //            msrDensityStep(msr,dTime);
    //        }
    //        
    //        if (msr->param.bDeltaAccelStep) {
    //            //FIXME Change logging from file to stderr
    //            //fprintf(stdout,"Adding DeltaAccelStep dt\n");
    //            if (!msr->param.bDeltaAccelStepGasTree) {
    //                msrActiveType(msr,TYPE_ALL,TYPE_TREEACTIVE);
    //                msrBuildTree(msr,0,-1.0,1);
    //            }
    //            else {
    //                msrActiveType(msr,TYPE_GAS,TYPE_TREEACTIVE);
    //                msrBuildTree(msr,0,-1.0,1);
    //            }
    //            msrActiveType(msr,TYPE_ALL,TYPE_TREEACTIVE|TYPE_SMOOTHACTIVE);
    //            msrBuildTree(msr,0,-1.0,1);
    //            msrActiveType(msr,TYPE_ALL,TYPE_ACTIVE|TYPE_TREEACTIVE|TYPE_SMOOTHACTIVE);
    //            /* This smooth sets dt directly -- hardwired coefficient */
    //            msrSmooth(msr,dTime,SMX_DELTAACCEL,0);
    //        }
    //        msrReorder(msr);
    //        //nOutputList = 0;
    //        //OutputList[nOutputList++]=OUT_DT_ARRAY;
    //        if(msr->param.iMaxRung > 1
    //           && (msr->param.iRungForceCheck || msr->param.bDensityStep || msrDoGravity(msr))) {
    //            msrActiveType(msr,TYPE_ALL,TYPE_ACTIVE|TYPE_TREEACTIVE|TYPE_SMOOTHACTIVE);
    //            msrDtToRung(msr,0,msrDelta(msr),1);
    //            msrRungStats(msr);
    //            //OutputList[nOutputList++]=OUT_RUNG_ARRAY;
    //        }
    //        
    //        //msrWriteOutputs(msr, achFile, OutputList, nOutputList, dTime);
    //        
    //        if (msr->param.iRungForceCheck != 0) {
    //            if (msr->param.iRungForceCheck > 0)
    //                msrActiveMaskRung(msr,TYPE_ACTIVE,msr->param.iRungForceCheck,1);
    //            else {
    //                assert(msr->param.iRungForceCheck > 0);
    //                /* Set active randomly -- e.g. every second particle */
    //                /* msriOrderToRung(msr,0,msrDelta(msr),1); */
    //            }
    //            
    //            
    //            printf("Force check on %d particles active on rung %d (max %d)\n",msr->nActive,msr->param.iRungForceCheck,msr->iCurrMaxRung);
    //            if (msr->param.iRungForceCheck > msr->iCurrMaxRung) msr->param.iRungForceCheck = msr->iCurrMaxRung;
    //            
    //            if (msr->nActive) {
    //                int nActive;
    //                msrDomainDecomp(msr,msr->param.iRungForceCheck,1);
    //                msrInitAccel(msr);
    //                
    //                printf("Forces, Step:%f nActive %i\n",0.,msr->nActive);
    //                if(msrDoGravity(msr)) {
    //                    if (msr->param.bDoSelfGravity) {
    //                        msrActiveRung(msr,msr->param.iRungForceCheck,1);
    //                        msrUpdateSoft(msr,dTime);
    //                        msrActiveType(msr,TYPE_ALL,TYPE_TREEACTIVE);
    //                        printf("Gravity, iRung: %d to %d\n", msr->param.iRungForceCheck, msr->param.iRungForceCheck);
    //                        msrBuildTree(msr,0,dMass,0);
    //                    }
    //                    msrGravity(msr,0.0,msrDoSun(msr),&iSec,&dWMax,&dIMax,&dEMax,&nActive);
    //                }
    //                
    //#ifdef GASOLINE
    //                if(msrDoGas(msr) && msrSphCurrRung(msr,msr->param.iRungForceCheck,1)) {
    //                    printf("SPH: iRung %d to %d\n",msr->param.iRungForceCheck,msr->param.iRungForceCheck);
    //                    msrSph(msr, dTime, msr->param.iRungForceCheck);
    //                }
    //#endif			 
    //                msrReorder(msr);
    //                //nOutputList = 2;
    //                //OutputList[0]=OUT_ACCELRFC_VECTOR;
    //                //OutputList[1]=OUT_PDVRFC_ARRAY;
    //                //if (msrDoDensity(msr)) OutputList[nOutputList++]=OUT_DENSITYRFC_ARRAY;
    //                
    //                //msrWriteOutputs(msr, achFile, OutputList, nOutputList, dTime);
    //            }
    //            
    //        }
    //    }
    
    //msr->param.iStartStep = msr->param.iStopStep;
    return 0;
}

int cleanup_code(){
    //dfFinalize( msr->df[0] );
    
    msrFinish(msr);
    mdlFinish(mdl);
    return 0;
}

int get_mass(int id, double * mass){
    return 0;
}

int commit_particles(){
    return 0;
}

int get_time(double * time){
    *time = dTime;
    return 0;
}

int set_mass(int id, double mass){
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
int new_sph_particle(int * index_of_the_particle, double mass, double x, 
                     double y, double z, double vx, double vy, double vz, 
                     double eps, double phi, double rho, double u,
                     double metals){
    msr->nGas += 1;
    msr->N += 1;
    msr->nMaxOrderGas = msr->nGas - 1;
    msr->nMaxOrder = NIORDERGASBUFFER + msr->N - 1;
    msr->nMaxOrderDark = NIORDERGASBUFFER + msr->nGas + msr->nDark - 1;
    
    *index_of_the_particle = msr->N;
    i = *index_of_the_particle;
    p = &pkd->pStore[i];
    
    //    p->iOrder = ???;
    //#if NIORDERGASBUFFER
    //    if (p->iOrder >= pkd->nGas) p->iOrder += NIORDERGASBUFFER;
    //#endif
    
    p->r[0] = x;
    p->r[1] = y;
    p->r[2] = z;
    p->v[0] = vx;
    p->v[1] = vy;
    p->v[2] = vz;
    p->fMass = mass;
    p->fSoft = eps;
    p->fPot = phi;
    p->fDensity = rho;
    p->u = u;
    p->uPred = u;
    p->fMetals = metals;
    p->fMetalsPred = metals;
    
    return 0;
}

int new_dm_particle(int * index_of_the_particle, double mass, double x, 
                    double y, double z, double vx, double vy, double vz, 
                    double eps, double phi){
    msr->nDark += 1;
    msr->N += 1;
    msr->nMaxOrderDark = NIORDERGASBUFFER + msr->nGas + msr->nDark - 1;
    msr->nMaxOrder = NIORDERGASBUFFER + msr->N - 1;
    
    p->r[0] = x;
    p->r[1] = y;
    p->r[2] = z;
    p->v[0] = vx;
    p->v[1] = vy;
    p->v[2] = vz;
    p->fMass = mass;
    p->fSoft = eps;
    p->fPot = phi;
    return 0;
}

int new_star_particle(int * index_of_the_particle, double mass, double x, 
                      double y, double z, double vx, double vy, double vz,
                      double eps, double phi, double metals, double tform){
    msr->nStar += 1;
    msr->N += 1;
    msr->nMaxOrder = NIORDERGASBUFFER + msr->N - 1;
    return 0;
}

int delete_particle(int index_of_the_particle){
    // Check if particle is dm, star or sph
    // make sure index_of_the_particle is adjusted for all following particles
    msr->N -= 1;
    return -1;
}

int get_potential(int index_of_the_particle, double * potential){
    return 0;
}

int synchronize_model(){
    return 0;
}

int set_state_sph(int index_of_the_particle, double mass, double x, double y, 
                  double z, double vx, double vy, double vz, double u){
    return 0;
}

int get_state_sph(int index_of_the_particle, double * mass, double * x, 
                  double * y, double * z, double * vx, double * vy, double * vz,
                  double * u){
    //id = msr->pMap[index_of_the_particle];
    return 0;
}

int set_state_star(int index_of_the_particle, double mass, double x, double y, 
                   double z, double vx, double vy, double vz){
    return 0;
}

int get_state_star(int index_of_the_particle, double * mass, double * x, 
                   double * y, double * z, double * vx, double * vy, double * vz
                   ){
    //id = msr->pMap[index_of_the_particle];
    return 0;
}

int set_state(int index_of_the_particle, double mass, double x, double y, 
              double z, double vx, double vy, double vz, double u){
    return 0;
}

int get_state(int index_of_the_particle, double * mass, double * x, 
              double * y, double * z, double * vx, double * vy, double * vz,
              double * u){
    //id = msr->pMap[index_of_the_particle];
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

int set_begin_time(double time){
    return 0;
}

int set_radius(int index_of_the_particle, double radius){
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

int get_velocity(int id, double * vx, double * vy, double * vz){
    return 0;
}

int get_position(int id, double * x, double * y, double * z){
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

int set_position(int id, double x, double y, double z){
    return 0;
}

int get_acceleration(int index_of_the_particle, double * ax, double * ay, 
                     double * az){
    return 0;
}

int commit_parameters(){
    return 0;
}

int set_velocity(int id, double vx, double vy, double vz){
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
