        SUBROUTINE GETPARTINFOC(NAME, LOCNUM, JTYPE, USER_NUMBER, JRCD)
        EXTERNAL GETPARTINFO
        write(*,*) "Fortran"
        CALL GETPARTINFO(LOCNUM, JTYPE, NAME, ELEM_NUMBER, JRCD)
        END