!  param_A23A.h
!
      integer(C_INT) npc,mx,my,mz,mxyz0,np0,myA,kd,   &
                     mxyz,mxyz3,mxyzA,                &
                     nob,nob2,nob3,iblk,iblk2,iblk3,  &
                     nhistm
!
      character :: praefixs*29,praefixc*23,suffix2*2,suffix1*2
!     character :: praefixs*29,praefixc*25,suffix2*2,suffix1*2
!
      parameter  (npc=8)
      parameter  (mx=72,my=144,mz=72)
      parameter  (mxyz0=mx*my*mz,np0=32*mx*my*mz) 
!
      parameter  (suffix2='0a',suffix1='0a') 
!     parameter  (suffix2='0b',suffix1='0a')
!                 ++++++++++++
!
      parameter (praefixs='/home/tanakam/mrg37/rec_3d23A',  &
                 praefixc='/data/sht/tanakam/forta')
!     parameter (praefixs='/home/mtanaka/mrg37/rec_3d23A',  &
!                praefixc='/home/mtanaka/mrg37/forta')
!
!            {ptx,pty,ptz} = (-2:mx+1,-2:my+2,-2:mz+1)
!                     {ex} = (-2:mx+1,-2:my+2,-2:mz+1)
!            myA: 0,1,...,my  start from 0
!
!  common/ptable/ gx(-2:mx+1),gy(-1:my+1),gz(-2:mz+1), 
      parameter  (myA=my+1,kd=mz/npc)           ! mz/npc must be divided
      parameter  (mxyz=mx*myA*mz,mxyz3=3*mx*myA*mz)
      parameter  (mxyzA=(mx+4)*(myA+2)*(mz+4))  ! -1,0,..,my,my+1
! 
      parameter  (nob=15,nob2=19,nob3=7) 
      parameter  (iblk=3,iblk2=1,iblk3=1)
      parameter  (nhistm=12)
