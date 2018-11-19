      subroutine errmess(lun11,nlen,str1) 
      character*(*) str1 
      integer lun11,nlen,ntmp
      ntmp=nlen
      write (lun11,*) 'in errmess:',str1 
      return 
      END                                           
