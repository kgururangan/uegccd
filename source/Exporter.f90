Module Exporter

Use Precision
Use Constants
Use Pointers

Implicit None

private
public :: ExportUEGIntegrals, ExportFCIDUMP

Contains
   
        Subroutine ConvertFockToHcore(HCore,HEGData,UEGInfo)
            
              Use Types, only: HEGDataType, UEGInfoType
              
              Type (HEGDataType), Intent(In) :: HEGData
              Type (UEGInfoType), Intent(In) :: UEGInfo
              Real (Kind=pr), Intent(Out) :: HCore(UEGInfo%NAO)
           
              Integer :: p, i
              Real (Kind=pr) :: g(UEGInfo%NAO), eri1, eri_x
              
              Do p=1,UEGInfo%NAO
                 !Hcore(p) = HEGData%Eigen(p)
                 g(p) = 0.0_pr
                 Do i=1,UEGInfo%NOcc
                    eri1 = ERI(HEGData,UEGInfo,p,i,p,i,DummyFlag=0)
                    eri_x = ERI(HEGData,UEGInfo,p,i,i,p,DummyFlag=0)
                    if (Min(eri1,eri_x) < -80_pr) cycle
                    g(p) = g(p) + 2.0_pr * eri1 - eri_x
                    !g(p) = g(p) - eri_x
                 End Do
                 HCore(p) = HEGData%Eigen(p) - g(p)
              End Do
              
        End Subroutine ConvertFockToHcore
   
        Subroutine ExportFCIDUMP(HEGData,UEGInfo)
              
              Use Types, only: HEGDataType, UEGInfoType
              Use HEG, only: FindIndex
              
              Type (HEGDataType), Intent(In) :: HEGData
              Type (UEGInfoType), Intent(In) :: UEGInfo

              Integer, Parameter :: io=101
              Integer :: p, q, r, s, cnt, a, b, i, j, x, y
              Real (Kind=pr) :: eri_val, eri_val2, fock_val
              Real (Kind=pr) :: HCore(UEGInfo%NAO)
              
              !Call ConvertFockToHcore(HCore,HEGData,UEGInfo)

              open(unit=io, file='FCIDUMP', status='replace')
              write(io,*) '&FCI ', "NORB=", UEGInfo%NAO, ",", "NELEC=", UEGInfo%Nelectron, ",", "MS2=", 0
              write(io,*) 'ORBSYM=', ','
              write(io,*) 'ISYM=1 UHF=.FALSE.'
              write(io,*) '&END'
              ! Write all 2e-integrals straight
              do p=1,UEGInfo%NAO
                 do q=1,UEGInfo%NAO
                    do r=1,UEGInfo%NAO
                       ! Use the momentum-allowed indexing
                       s = FindIndex(HEGData,UEGINfo%MaxKPoint,UEGInfo%NAO,p,q,r)
                       if (s<=0) cycle
                       eri_val = ERI(HEGData,UEGInfo,p,q,r,s,DummyFlag=0)
                       if (eri_val < -80_pr) cycle
                       write(io,*) eri_val, p, q, r, s
                    end do
                 end do
              end do
              ! Write only permutationally unique integrals
!              do p = 1,UEGInfo%NAO
!                 do r = 1,p
!                    x = p*(p + 1)/2 + r - 1
!                    do q = 1,UEGInfo%NAO
!                       !do s = 1,q
!                          s = FindIndex(HEGData,UEGINfo%MaxKPoint,UEGInfo%NAO,p,q,r)
!                          if (s<=0) cycle
!                          y = q*(q + 1)/2 + s - 1
!                          if (x >= y) then
!                              eri_val = ERI(HEGData,UEGInfo,p,q,r,s,DummyFlag=0)
!                              if (eri_val < -80_pr) cycle
!                              write(io,*) eri_val, p, r, q, s
!                          end if
!                       !end do
!                    end do
!                 end do
!              end do
              ! Write 1e-integrals
              do p = 1,UEGInfo%NAO
                 do q = 1,UEGInfo%NAO
                    if (p == q) then
                       !fock_val = HCore(p)
                       fock_val = HEGData%Eigen(p)
                    else
                       fock_val = 0.0_pr
                    end if
                    write(io,*) fock_val, p, q, 0, 0
                 end do
              end do
              ! write nuclear repulsion (set to Madelung constant)
              write(io,*) HEGData%Madelung, 0, 0, 0, 0
              close(io)
      End Subroutine ExportFCIDUMP

      Subroutine ExportUEGIntegrals(HEGData,UEGInfo)
              
              Use Types, only: HEGDataType, UEGInfoType
              Use HEG, only: FindIndex
              
              Type (HEGDataType), Intent(In) :: HEGData
              Type (UEGInfoType), Intent(In) :: UEGInfo

              Integer, Parameter :: io=101
              Integer :: p, q, r, s, cnt, a, b, i, j
              Real (Kind=pr) :: eri_val, eri_val2, fock_val
              Real (Kind=pr) :: HCore(UEGInfo%NAO)
              
              !Call ConvertFockToHcore(HCore,HEGData,UEGInfo)

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
