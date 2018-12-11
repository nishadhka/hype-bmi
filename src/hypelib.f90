MODULE HYPELIB

    USE MODELMODULE, ONLY :         model,                                &
                                    model_version_information,            &
                                    define_output_variables,              &
                                    set_modelconfig_from_parameters,      &
                                    initiate_model,                       &
                                    initiate_model_state,                 &
                                    define_model_parameters,              &
                                    load_modeldefined_input
    USE STATETYPE_MODULE
    USE WORLDVAR, ONLY :            writematlab,                          &
                                    writeload,                            &
                                    output,                               &
                                    noutput,                              &
                                    nacrit,                               &
                                    nsubCrit,                             &
                                    ndt,                                  &
                                    fileunit_temp,                        &
                                    fileunit_tests,                       &
                                    maxcharpath,                          &
                                    simsequence,                          &
                                    infodir,                              &
                                    modeldir,                             &
                                    resdir,                               &
                                    forcingdir,                           &
                                    logdir,                               &
                                    bdate,                                &
                                    sdate,                                &
                                    outstartdate,                         &
                                    dtskip,                               &
                                    doopt,                                &
                                    numoptimpar,                          &   
                                    deallocate_worldvar,                  &
                                    deallocate_MCvariables,               &
                                    optim,                                &
                                    bestMCoptcrit,                        &
                                    bestMCperformance,                    &
                                    bestMCparameters,                     &
                                    maxperf,                              &
                                    maxsubass,                            &
                                    simsubmodel,                          &
                                    ibasemodel,                           &
                                    psdates,                              &
                                    optimStartTime,                       &      
                                    optimFuncCall,                        &
                                    lineSearchCallCount,                  &
                                    doassimilation,                       &
                                    noutreg,                              &
                                    da_allocate_accumulation,             &
                                    allocate_accumulation,                &
                                    reallocate_outvar_information,        &
                                    usestop84,                            &
                                    indatacheckonoff,                     &
                                    indatachecklevel
    USE MODVAR, ONLY :              nsub,                                 &
                                    ncrop,                                &
                                    nsub_basemodel,                       &
                                    nclass,                               &
                                    numsubstances,                        &
                                    naquifers,                            &
                                    maxsoillayers,                        &
                                    max_classoutvar,                      &
                                    max_basinoutvar,                      &
                                    max_noutvar,                          &
                                    conduct,                              &
                                    statesize,                            &
                                    conductxoms,                          &
                                    conductregest,                        &
                                    conductwarning,                       &
                                    simulate,                             &
                                    allocate_outvar,                      &
                                    deallocate_modvar,                    &
                                    currentdate,                          &
                                    dayno,                                &
                                    noutvar,                              &
                                    noutvarclass,                         &
                                    nrivertypes,                          &
                                    nlaketypes,                           &
                                    outvar
    USE COMPOUT, ONLY :             compute_mapoutput,                    &
                                    compute_outloads,                     &
                                    prepare_to_compute_crit,              &
                                    calculate_criteria
    USE TIMEROUTINES, ONLY :        calculate_time_for_model
    USE READWRITE_ROUTINES
    USE LIBDATE, ONLY :             DateType, OPERATOR(.EQ.)
    USE DATAMODULE, ONLY:           get_hyss_arguments
    USE OPTIMIZATION
    USE STATE_DATAMODULE, ONLY :    initiate_state_for_submodel,          &
                                    load_saved_state,                     &
                                    finalize_outstate
    USE MODEL_TEST_ROUTINES, ONLY : run_hype_setup_tests,                 &
                                    run_hype_finalize_tests,              &
                                    run_hype_tests,                       &
                                    run_hype_observation_tests

    IMPLICIT NONE
    ! Parameter declarations
    INTEGER, PARAMETER :: maxoutstates = 10     !Max number of dates for saving state

    ! Variable declarations
    TYPE(DateType)                  d                           !Current time
    INTEGER                         idt,        &               !Current timestep
                                    ivar,       &               !Current output variable
                                    iens,       &               !Current ensemble being simulated
                                    iout                        !Current output
!   INTEGER nobutton                                            !No exit window
    LOGICAL ::                      pwrite                      !Flag for periodend, time to write to file
    CHARACTER(LEN=maxcharpath+25)   filename                    !hyss filename
    CHARACTER(LEN=8)  ::            logdate                     !Date for log-file name
    CHARACTER(LEN=10) ::            logtime                     !Time for log-file name
    CHARACTER(LEN=3)  ::            logseq                      !Seqnr for log-file name
    INTEGER ::                      datim(8),   &               !Date in YYYYmmdd
                                    oldyear                     !year of last time step

    REAL, ALLOCATABLE ::            par(:)
    REAL                            optcrit,    & 
                                    condcrit,   &
                                    condthres
    REAL, ALLOCATABLE ::            basincrit(:,:,:)            !R2, CC, RE, RSDE, QC, QR, STDC, STDR, MAE, RMSE, Bias, STDbias, KGE, KGEpartSTD, KGEpartMM, NRMSE per subbasin och kriterie
    REAL, ALLOCATABLE ::            simperformance(:,:)         !rr2,sr2,wr2,rmae,sbias,rrve,wrve,rra,sra,meanRA,tau,medianr2,medianra,meanrs,meancc,mediankg,meanabsre
    INTEGER, ALLOCATABLE ::         subincrit(:)                !Subbasins to be included in criteria calculations
    INTEGER ::                      npar,       &               !Number of parameters to be calibrated (couted in file)
                                    nmapperiod, &               !Number of periods for map print out
                                    numoutstates                !Number of dates for saving state
    LOGICAL ::                      stateinput                  !Code for reading state
    TYPE(DateType) ::               stateoutadate(maxoutstates)       !Dates for saving state
    TYPE(DateType) ::               prestateoutdate(maxoutstates)     !Day before date for saving state (or 0)
      
    ! Variables for updating of Q and W
    LOGICAL ::                      quseobsallstations,         &
                                    quseobsnostations,          &
                                    qarnostations,              & 
                                    warnostations,              &
                                    wendupdallstations,         & 
                                    wendupdnostations
    CHARACTER(LEN=4) ::             wobsvarname
      
    ! Model state variables and other saved variables declaration
    TYPE(SNOWICESTATETYPE) ::       frozenstate
    TYPE(SOILSTATETYPE)    ::       soilstate      
    TYPE(AQUIFERSTATETYPE) ::       aquiferstate      
    TYPE(RIVERSTATETYPE)   ::       riverstate      
    TYPE(LAKESTATETYPE)    ::       lakestate
    TYPE(MISCSTATETYPE)    ::       miscstate

CONTAINS

    FUNCTION initialize(dir, iseq) RESULT(istat)

        CHARACTER(LEN=maxcharpath), INTENT(IN), OPTIONAL    :: dir
        INTEGER, INTENT(IN), OPTIONAL                       :: iseq
        INTEGER                                             :: istat
        
        CALL DATE_AND_TIME(logdate, logtime, values=datim)
        
        CALL model_version_information(0)

        IF(PRESENT(dir)) THEN
            infodir=dir
        ELSE
            CALL getcwd(infodir, istat)
        ENDIF
        IF(PRESENT(iseq)) THEN
            simsequence=iseq
        ENDIF

        WRITE(logseq, '(I3.3)') simsequence
        WRITE(filename, '(a)') TRIM(infodir)//'hyss_'//logseq(1:3)//'_'//logdate(3:8)//'_'//logtime(1:4)//'.log'
        OPEN(UNIT=6,FILE=TRIM(filename),STATUS = 'replace',ACTION='write')
        WRITE(6, '(A,I4,A,I2.2,A,I2.2,A,I2.2,A,I2.2,A,I2.2)')              &
                 'Job start date: ',datim(1),'-',datim(2),'-',datim(3),    &
                 'time: ',datim(5),':',datim(6),':',datim(7)
        WRITE(6, *) '---------------------------------------------------'
!       nobutton = SETEXITQQ(qwin$exitnopersist)    !nobutton version
        CALL model_version_information(6)   !Print model version information in logg-file
        
        CALL define_output_variables()
        
        CALL define_model_parameters()
        
        CALL load_coded_info(infodir, istat, bdate, sdate, outstartdate, dtskip, ndt, numsubstances,  &
                             stateinput, maxoutstates, numoutstates, stateoutadate, prestateoutdate,   &
                             modeldir, resdir, forcingdir, logdir,                                     &
                             quseobsallstations, quseobsnostations, qarnostations, warnostations,      &
                             wendupdallstations, wendupdnostations, wobsvarname, subincrit)
      IF(istat/=0) RETURN
      nmapperiod = 0 !initialization needed if no mapoutput

      CALL load_basindata(modeldir, forcingdir, nsub_basemodel, istat)
      IF(istat.NE.0) RETURN

      CALL load_cropdata(modeldir, ncrop, istat)
      IF(istat.NE.0) RETURN

      CALL load_pointsourcedata(modeldir, 'PointSourceData.txt', nsub_basemodel, istat) 
      IF(istat.NE.0) RETURN
      
      CALL load_soilleakage_concentrations(modeldir, nsub_basemodel, istat)
      IF(istat.NE.0) RETURN
      
      CALL load_branchdata(modeldir, istat)
      IF(istat.NE.0) RETURN
      
      CALL load_aquiferdata(modeldir, nsub_basemodel, naquifers, istat) 
      IF(istat.NE.0) RETURN

      CALL load_glacierdata(modeldir, nsub_basemodel, istat) 
      IF(istat.NE.0) RETURN

      CALL initiate_model_parameters(nsub_basemodel, istat)
      IF(istat.NE.0) RETURN

      WRITE(filename,'(a)') TRIM(logdir)//'tests_'//logseq(1:3)//'_'//logdate(3:8)//'_'//logtime(1:4)//'.log'
      CALL run_hype_setup_tests(fname=filename, onoff=indatacheckonoff, level=indatachecklevel) !testetup

      CALL run_hype_observation_tests(istat)
      IF(istat.NE.0) RETURN

      CALL load_submodel_info(infodir, simsubmodel, nsub, istat)   !allocation and initialisation of ibasemodel
      IF(istat.NE.0) RETURN

      CALL load_output_regions(modeldir, noutreg, istat)    !Prepare for outregions, read Outregions.txt
      IF(istat.NE.0) RETURN
      
      CALL load_observations(forcingdir, nsub_basemodel, bdate, sdate, ndt, istat)
      IF(istat.NE.0) RETURN
      
      CALL load_parameters(modeldir, nsub_basemodel, 'par.txt')
      
!Read base model defined input data and set base model configuration
      CALL set_model_base_configuration(nsub_basemodel, stateinput, modeldir, forcingdir,&
                                      & statesize, istat)
      IF(istat.NE.0) RETURN
      
      IF(simsubmodel) THEN
        CALL reform_inputdata_for_submodel(nsub_basemodel, nsub, ibasemodel)
      ENDIF
      CALL calculate_path(nsub)

!Read model defined input data and set model configuration
      CALL load_modeldefined_input(modeldir, forcingdir, nsub_basemodel, nsub, ibasemodel,&
                                 & bdate, sdate, conductxoms, conductregest, istat)
      IF(istat.NE.0) RETURN

      CALL set_model_configuration(conduct, simulate)

      CALL set_modelconfig_from_parameters()
      
!Test input data and validate options
      CALL run_hype_tests(istat)

      CALL run_hype_finalize_tests()
      
      IF(istat.NE.0) RETURN

      CALL prepare_for_update(modeldir, wobsvarname, quseobsallstations, quseobsnostations,
                            & qarnostations,warnostations,wendupdallstations,wendupdnostations,
                            & nsub)

!>Initialisations for memory allocation (states and output)
      CALL allocate_model_states(nsub_basemodel,statesize,conduct, frozenstate, soilstate, &
                               & aquiferstate, riverstate, lakestate, miscstate)

      !CALL set_output_default_classgroup(noutput)
      istat = set_outvar(noutput, noutvar, noutvarclass)
      IF(istat.NE.0) RETURN

      istat = set_outvar_crit(noutvar, noutvarclass)
      IF(istat.NE.0) RETURN

      IF(.NOT.doassimilation) THEN
        CALL allocate_outvar(noutvar, noutvarclass)
        CALL reallocate_outvar_information(noutvar, noutvarclass)
        CALL allocate_accumulation(nsub, nclass, numsubstances, max_classoutvar, max_basinoutvar)
      ENDIF

!Allocate local variables
      nsubCrit = nsub
      ALLOCATE(basincrit(nsubCrit, maxsubass, nacrit))
      ALLOCATE(simperformance(maxperf, nacrit))

!Preparations for subbasin output
      CALL prepare_subbasin_output(subincrit, istat)
      IF(istat.NE.0) RETURN
      
      CALL DATE_AND_TIME(values=datim)
      WRITE(6,*)
      WRITE(6,*) '---------------------------------------------------'
      WRITE(6,'(A,I4,A,I2.2,A,I2.2,A,I2.2,A,I2.2,A,I2.2)')    &
          ' Initialisations finished, calculations starts: ',datim(1),'-',datim(2),'-',datim(3),    &
          '  time: ',datim(5),':',datim(6),':',datim(7)
      WRITE(6,*) '---------------------------------------------------'
      FLUSH 6

      !>Simulation start
      !>Initial model calculations; initial states, parameters
      CALL prepare_outputfiles(resdir, nsub, naquifers, 1, .FALSE., .FALSE.)
      CALL initiate_output_routines()     !All output accumulation variables zeroed
      CALL set_modelconfig_from_parameters()
      IF(simsubmodel)THEN
          CALL initiate_state_for_submodel(forcingdir, ibasemodel, stateinput, frozenstate, soilstate,&
              & aquiferstate, riverstate, lakestate, miscstate)
      ELSE
          IF(stateinput) THEN
              WRITE(6,*) 'Loading saved state.'
              CALL load_saved_state(forcingdir, nsub, frozenstate, soilstate, aquiferstate, riverstate, &
                  & lakestate, miscstate)
          ELSE
              CALL initiate_model_state(frozenstate, soilstate, aquiferstate, riverstate, lakestate, miscstate)
          ENDIF
          CALL initiate_model(frozenstate, soilstate, aquiferstate, riverstate, lakestate, miscstate) 
      ENDIF
      oldyear = 0

    END SUBROUTINE initialize


    END MODULE HYPELIB
