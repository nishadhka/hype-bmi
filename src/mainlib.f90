PROGRAM MAINLIB

    USE HYPELIB, ONLY:      initialize, update, finalize, get_num_subbasins, get_num_output_fields,&
                            get_latlons, get_output_fields, get_basin_field
    USE WORLDVAR, ONLY:     maxcharpath, ndt
    USE DATAMODULE, ONLY:   get_hyss_arguments

    IMPLICIT NONE

    CHARACTER(LEN=maxcharpath)                          :: dir
    INTEGER                                             :: istat, iseq, idt, count, ifld, nflds
    CHARACTER(64), ALLOCATABLE                          :: fldids(:)
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
    CALL get_num_output_fields(nflds)
    ALLOCATE(fldids(nflds))
    CALL get_output_fields(fldids)
    DO idt = 1, ndt
        istat = update()
        IF(istat/=0) STOP istat
        WRITE(0,*) "------------- TS ", idt, "-------------"
        DO ifld = 1, nflds
            CALL get_basin_field(Q, fldids(ifld))
            WRITE(0,*) fldids(ifld), "max = ", MAXVAL(Q), ", min = ", MINVAL(Q)
        END DO
    ENDDO
    istat = finalize()
    IF(istat/=0) STOP istat
    DEALLOCATE(lats, lons, Q)
    STOP 0

END PROGRAM
