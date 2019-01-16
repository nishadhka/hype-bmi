PROGRAM MAINLIB

    USE HYPELIB, ONLY:      initialize, update, finalize, get_num_subbasins, get_latlons, get_discharge
    USE WORLDVAR, ONLY:     maxcharpath, ndt
    USE DATAMODULE, ONLY:   get_hyss_arguments

    IMPLICIT NONE

    CHARACTER(LEN=maxcharpath)                          :: dir
    INTEGER                                             :: istat, iseq, idt, count
    REAL, ALLOCATABLE                                   :: lats(:), lons(:), Q(:)

    CALL get_hyss_arguments(dir, iseq)
    istat = initialize(dir, iseq)
    count = get_num_subbasins()
    WRITE(0,*) "The number of subbasins is", count
    ALLOCATE(lats(count), lons(count), Q(count))
    IF(istat/=0) STOP istat
    CALL get_latlons(lats, lons)
    WRITE(0,*) "The latitudes are",  lats
    WRITE(0,*) "The longitudes are", lons
    DO idt = 1, ndt
        istat = update()
        IF(istat/=0) STOP istat
        CALL get_discharge(Q)
        WRITE(0,*) "------------- TS ", idt, "-------------"
        WRITE(0,*) "Q-max = ", MAXVAL(Q)
    ENDDO
    istat = finalize()
    IF(istat/=0) STOP istat
    DEALLOCATE(lats, lons, Q)
    STOP 0

END PROGRAM
