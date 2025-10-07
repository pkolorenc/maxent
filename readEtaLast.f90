PROGRAM ReadEta

  INTEGER iou,NT,maxNT
  LOGICAL addLog
  REAL(8), ALLOCATABLE :: eta(:)
  CHARACTER(4) :: arg

  call get_command_argument(1,arg)
  read(arg,*) maxNT

  iou = 1000
  open(iou,file='eta.last',status='old', &
    access='sequential',form='unformatted')
  read(iou) NT,addLog
  if (addLog) then
    allocate(eta(0:maxNT+1))
  else
    allocate(eta(0:maxNT))
  end if
  eta(:) = 0.d0
  read(iou) eta(0:NT)
  close(iou)

  if (addLog) then
    eta(maxNT+1) = eta(NT)
    eta(NT) = 0.d0
  end if

  write(6,'(I4," ")',advance='no') NT
  do iou = 0,size(eta)-1
    write(6,'(" ",ES20.10)',advance='no') eta(iou)
  end do
  write(6,*)

  deallocate(eta)

END PROGRAM
