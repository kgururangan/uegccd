Module Exporter

Use Precision
Use Constants
Use Pointers

Implicit None

private
public :: ExportUEGIntegrals

Contains

        Subroutine ExportUEGIntegrals(HEGData,UEGInfo)
              
              Use Types, only: HEGDataType, UEGInfoType
              Use HEG, only: FindIndex
              
              Type (HEGDataType), Intent(In) :: HEGData
              Type (UEGInfoType), Intent(In) :: UEGInfo

              Integer, Parameter :: io=101
              Integer :: p, q, r, s, cnt, a, b, i, j
              Real (Kind=pr) :: eri_val, eri_val2, fock_val

              write(6,*) 'Printing UEG Info Parameters'
              write(6,*) '---------------------------'
              write(6,*) 'Number of Electrons:', UEGInfo%Nelectron
              write(6,*) 'Number of Occupied:', UEGInfo%NOcc
              write(6,*) 'Number of basis functions:', UEGInfo%NAO

              open(unit=io, file='ueg.inp', status='replace')
              write(io,*) UEGInfo%Nelectron, UEGInfo%NAO, UEGInfo%NOcc
              write(io,*) HEGData%EHF
              close(io)

              open(unit=io, file='onebody.inp', status='replace')
              write(6,*) 'Writing onebody.inp'
              cnt = 1
              do p = 1,UEGInfo%NAO
                 do q = 1,p
                    if (p == q) then
                       fock_val = HEGData%Eigen(p)
                    else
                       fock_val = 0.0_pr
                    end if
                    write(io,*) fock_val, cnt
                    cnt = cnt + 1
                 end do
              end do
              close(io)
              open(unit=io, file='twobody.inp', status='replace')
              write(6,*) 'Writing twobody.inp'
              ! Write only permutationally unique integrals
              !do p = 1,nbf
              !   do r = 1,p
              !      pr = p*(p + 1)/2 + r
              !      do q = 1,nbf
              !         do s = 1,q
              !            qs = q*(q + 1)/2 + s
              !            if (pr > kl) then
              !                write(io,*) p, q, r, s, eri_val
              !            end if
              !         end do
              !      end do
              !   end do
              !end do
              ! Write all integrals straight
              do p=1,UEGInfo%NAO
                 do q=1,UEGInfo%NAO
                    do r=1,UEGInfo%NAO
                       ! Use the momentum-allowed indexing
                       s = FindIndex(HEGData,UEGINfo%MaxKPoint,UEGInfo%NAO,p,q,r)
                       if (s<=0) cycle
                       eri_val = ERI(HEGData,UEGInfo,p,q,r,s,DummyFlag=0)
                       if (eri_val < -80_pr) cycle
                       write(io,*) p, q, r, s, eri_val
                    end do
                 end do
              end do
              ! write nuclear repulsion (set to Madelung constant)
              write(io,*) 0, 0, 0, 0, HEGData%Madelung
              close(io)
              !open(unit=io, file='mask', status='replace')
              !write(6,*) 'Writing T2 mask'
              !do a=nocc+1,nbf
              !   do i=1,nocc
              !      do j=1,nocc
              !         b = FindIndex(i,j,a)
              !         if (b <= nocc) cycle
              !         write(io,*) a - nocc, b - nocc, i, j
              !      end do
              !   end do
              !end do
              !close(io)
      End Subroutine ExportUEGIntegrals

End Module Exporter
