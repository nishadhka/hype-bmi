PROGRAM MAINLIB

    USE HYPELIB, ONLY:      initialize, update, finalize
    USE WORLDVAR, ONLY:     maxcharpath, ndt
    USE DATAMODULE, ONLY:   get_hyss_arguments

    IMPLICIT NONE

    CHARACTER(LEN=maxcharpath)                          :: dir
    INTEGER                                             :: istat, iseq, idt

    CALL get_hyss_arguments(dir, iseq)
    istat = initialize(dir, iseq)
    IF(istat.NE.0) STOP istat
    DO idt = 1, ndt
        istat = update()
        IF(istat.NE.0) STOP istat
    ENDDO
    istat = finalize()
    IF(istat.NE.0) STOP istat
    STOP 0

END PROGRAM
