module io
implicit none

character(len=:),allocatable :: infile, tFormat, nFormat
integer :: inunit

end module io


module simulation
implicit none
integer :: nTraj, iTraj   !number of trajectories and current trajectory
integer :: iSave          !the index of the most recent saved state (equal dt * int(t/dt) )
integer :: dtRecip,nSave  !number of steps between saved states, number of states to save
integer*8 :: nCritical      !threshold population for classifying a reaction as critical (more complex than this if non-first-order)
real*8  :: tSim, dt, t    !length of each trajectory, time interval between saves, current time, current stage
real*8  :: tauLeap, tauCritical !tauLeap is the actual step taken. tauCritical is the time until the next critical rxn fires
real*8  :: tauTentative   !tentative tau leap value
real*8  :: eps            !maximum tolerated estimated change in reaction propensities during tau leap. 0.03d0 is typical
real*8  :: propensitySum  !total reaction propensity for the whole system
real*8  :: tmpPropensitySum !for choosing which reaction fires
real*8  :: criticalPropensitySum !total reaction propensity for critical reactions
real*8  :: thresholdSSA   !determines when we switch to full SSA. Lower is more cautious. 10 is typical.
real*8  :: dtTolerance    !very small number that defines tolerance for digital precision loss (on the order of 10^-15)
real*8,allocatable,dimension(:) :: rates, mu, sigmaSq, maxEpsN !the current array of reaction rate constants, mean and variance of rates for each species
real*8,allocatable,dimension(:) :: propensities !the current array of reaction propensities
integer*8,allocatable,dimension(:,:) :: currentTrajectory !store the populations over time for each species
integer*8,allocatable,dimension(:) :: l !estimated number of times each reaction can fire without exhausting one of its reactants
!character(len=40) :: infile, tFormat, nFormat
!integer :: inunit
integer :: criticalReaction !the next reaction to fire without tau leaping
integer :: nSSAstep       !number of steps to run when reverting to full SSA, before checking for noncritical reactions again
logical,allocatable,dimension(:) :: isCritical !list of critical reactions
logical :: fullSSA,fireCritical,updateStage !true if tau is larger than the (random) time to the next critical rxn. if false no critical rxns fire
end module simulation



module model
implicit none
integer :: nSpecies, nRxn, nStages, stage                !number of species, number of reactions, number of stages (how many times to the rates get set?)
integer*8,allocatable,dimension(:)   :: n           !current populations
integer*8,allocatable,dimension(:)   :: reactants   !identifies which species is the reactant in each reaction (will have to modify
                                                  !this at some point, if I want to include more complex reaction types (non-first
                                                  !order)
integer,allocatable,dimension(:,:) :: v           !state change vectors for reactions
integer*8,allocatable,dimension(:,:) :: n0        !initial populations
real*8,allocatable,dimension(:)    :: switchTimes !when each stage *ends* last switch time must be larger than tSim
real*8,allocatable,dimension(:,:)  :: rateArray   !array of rates (each row contains the rates for one reaction, for all stages)
end module model








program stochastic
use simulation
use model
use io
implicit none

!function to calculate candidate tau leap
real*8 tentativeLeap, pro
!dummy array to assist in reading input file
real*8,allocatable,dimension(:) :: tmpRates
real*8,dimension(2) :: rand
!indices for looping over species and reactions
integer :: s, rxn
!function to draw from a Poisson distribution
integer*8 :: poidev
!generic iterator
integer :: i, remainingSSAsteps
!input file

infile = '/home/allanfriesen/local/stochastic/stochastic.in'
inunit = 10

thresholdSSA = 10.d0
nCritical = 10
dtTolerance = 1.d-12
nFormat = '(I22)'

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!                                      !!
!!  READ INPUT FILE AND ALLOCATE ARRAYS !!
!!                                      !!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


open(unit=inunit,file=infile)

!number of trajectories to generate
read(inunit,*) nTraj
read(inunit,*) tSim
read(inunit,*) dtRecip


!! inputs defining model and initial conditions !!

!number of unique types of "particles"
read(inunit,*) nSpecies

!arrays to hold populations and initial populations
allocate( n0(nTraj, nSpecies) )
allocate( n(nSpecies), mu(nSpecies), sigmaSq(nSpecies), maxEpsN(nSpecies) )

!read initial populations
do iTraj=1, nTraj
    read(inunit,*) (n0(iTraj,s), s=1, nSpecies)
end do !s=1, nSpecies

!number of reactions and stages (how many time intervals to define seperate rates)
read(inunit,*) nRxn
read(inunit,*) nStages

!arrays to hold state change vectors, rates, and switching times (when to change rates)
allocate( v(nRxn,nSpecies) )
allocate( rates( nRxn ), reactants(nRxn), l(nRxn), isCritical(nRxn), propensities(nRxn) )
allocate( rateArray(nRxn, nStages) )
allocate( switchTimes(nStages), tmpRates(nStages) ) !the switch time for the final stage just needs to be larger than tSim
!read state change vectors
do rxn=1, nRxn
    do s=1, nSpecies
        read(inunit,*) v(rxn,s)
    end do !s=1, nSpecies
    read(inunit,*) reactants(rxn)
end do !rx=1, nRxn
!read rates
do rxn=1, nRxn
    read(inunit,*) tmpRates
    do stage=1, nStages
        rateArray(rxn,stage) = tmpRates(stage)
    end do !i=1, nStages
end do! rx=1, nRxn
!read switch times
read(inunit,*) switchTimes

close(inunit)


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!                           !!
!! FINISH READING INPUT FILE !!
!!                           !!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!initialize some items needed for doing the simulation
nSave = tSim * dtRecip ! number of states to record in trajectory file *after* the initial state
dt = tSim / real(nSave,kind=8)
nSSAstep = 100 !number of SSA steps to execute when all reactions are critical
eps = 0.03d0  !maximum tolerated average relative change in rxn propensity during a tau leap. 0.03d0 is typical.
tFormat = '(F10.2)'
allocate( currentTrajectory(nSpecies, 0:nSave) )
!first row of output is the times when the system's state will be saved


do i=0, nSave
  write(*,tFormat,advance='no') i * dt
end do !i=0, nSave
write(*,*) !move to the next record in the trajectory file

!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!                        !!
!! LOOP OVER TRAJECTORIES !!
!!                        !!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                          !!
do iTraj=1, nTraj !!!!!!!!!!
                          !!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!



currentTrajectory(:,:) = 0.d0
fullSSA = .false.
do s=1, nSpecies
    currentTrajectory(s,0) = n0(iTraj,s)
    n(s) = n0(iTraj,s)
end do !s=1, nSpecies
do rxn=1, nRxn
    rates(rxn) = rateArray(rxn,1)
end do !rxn=1, nRxn
stage = 1
iSave = 0
t = 0.d0
remainingSSAsteps = nSSAstep



do while (t < tSim)
    updateStage = .false.
    call calcPropensities
    tauCritical = 0.d0
    if( remainingSSAsteps > 0 ) then
        call random_number(rand)
        tauCritical = log(1.d0 / rand(1))/propensitySum
        if( t + tauCritical > switchTimes(stage) ) then
            tauCritical = switchTimes(stage) - t
            t = switchTimes(stage)
            stage = stage + 1
            do rxn=1, nRxn
                rates(rxn) = rateArray(rxn,stage)
            end do !rxn=1, nRxn
            if( int(t / dt) > iSave) then
                call saveState !Need to fix this and/or saveState routine, so that state is saved through t + tauCritical
                               !For example, perhaps I could send the step as an argument...wait - not this one (see below)
            end if !int(t / dt) > iSave
            cycle
        end if !t + tau > switchTimes(state)
        t = t + tauCritical
        if (int(t/dt + dtTolerance) > iSave) then !modified this call on 13 Nov 2024, to avoid rare double calls
        !if( t > (iSave + 1) * dt) then
            call saveState
        end if !t + tau > (iSave + 1) * dt
        !t = t + tauCritical !THIS ONE -- I BELIEVE THE TIME SHOULD BE UPDATED BEFORE THE SAVESTATECALL
        tmpPropensitySum = 0.d0
        rand(2) = rand(2) * propensitySum
        criticalReaction = -1
        do rxn = 1, nRxn
            tmpPropensitySum = tmpPropensitySum + propensities(rxn)
            if( tmpPropensitySum > rand(2) ) then
                criticalReaction = rxn
                exit
            end if
        end do !
        if( criticalReaction < 0 ) then
            STOP 'Something is wrong: failed to select critical reaction during full SSA.' 
        end if !checks that you did indeed find a reaction to fire

        !update the species involved in the critical reaction
        
        do s = 1, nSpecies
            n(s) = n(s) + v(criticalReaction, s)
        end do !s = 1, nSpecies
        remainingSSAsteps = remainingSSAsteps - 1
        cycle
    else
        tauTentative = tentativeLeap()
        if( tauTentative + t > (iSave + 1) * dt ) then !don't leap past timestep for saving states
            tauTentative = (iSave + 1) * dt - t + dtTolerance
        end if ! tauTentative + t > (iSave + 1) * dt
        if( tauTentative + t > switchTimes(stage) ) then
            tauTentative = switchTimes(stage) - t + dtTolerance
            updateStage = .true.
        end if
        if( tauTentative > thresholdSSA / propensitySum ) then !do full SSA if at or *below* SSA
            call random_number(rand)  !declare rand --> probably will just need to be real*8(2), but make sure
            tauCritical = log(1.d0 / rand(1))/criticalPropensitySum !declare tauCritical
            if( tauTentative < tauCritical ) then
                tauLeap = tauTentative
                !no critical reactions will fire - just do the leap for non-critical reactions
            else
                if( tauTentative > tauCritical ) then
                    updateStage = .false.
                end if !tauTentative > tauCritical
                tauLeap = tauCritical
                !Use the random number to identify the critical reaction that fires
                tmpPropensitySum = 0.d0
                rand(2) = rand(2) * criticalPropensitySum
                criticalReaction = -1
                do rxn = 1, nRxn
                    if( isCritical(rxn) ) then
                        tmpPropensitySum = tmpPropensitySum + propensities(rxn)
                        if( tmpPropensitySum > rand(2) ) then
                            criticalReaction = rxn
                            exit
                        end if
                    end if
                end do !
                if( criticalReaction < 0 ) then
                    STOP 'Something is wrong: failed to select critical reaction.' 
                end if !checks that you did indeed find a reaction to fire
                !update the species involved in the critical reaction
                do s = 1, nSpecies
                    n(s) = n(s) + v(criticalReaction, s)
                end do !s = 1, nSpecies
                !update non-critical reactions with tau leaping
            end if !tauTentative < tauCritical
            do rxn = 1, nRxn
                if( .not.(isCritical(rxn)) ) then
                    pro = propensities( rxn )
                    do s = 1, nSpecies
                        n(s) = n(s) + poidev( pro * tauLeap ) * v(rxn,s)
                    end do !s = 1, nSpecies
                end if !.not. isCritical(rxn)
            end do !rxn = 1, nRxn
            t = t + tauLeap
            if( updateStage ) then
                t = switchTimes(stage)
                stage = stage + 1
                do rxn=1, nRxn
                    rates(rxn) = rateArray(rxn,stage)
                end do !rxn=1, nRxn
                updateStage = .false.
            end if !updateStage
            if( int(t/dt+dtTolerance) > iSave ) then
                call saveState
            end if !int(t/dt) > iSave

        else
            remainingSSAsteps = nSSAstep
        end if !tauTentative > thresholdSSA > propensitySum
    end if !remainingSSAsteps > 0
 
end do !while t < tSim




call saveState

do s=1, nSpecies
    do iSave=0,nSave
        write(*,nFormat,advance='no') currentTrajectory(s, iSave)
    end do !iSave=1, nSave
    write(*,*)
end do !s=1, nSpecies





!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
                              !!
end do !iTraj=1, nTraj !!!!!!!!!
                              !!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!                            !!
!! END LOOP OVER TRAJECTORIES !!
!!                            !!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!




end program





subroutine saveState
use model
use simulation
implicit none
integer :: s, countSave !countSave = DEBUGGING
do while (iSave < int(t/dt))
    countSave = 0 !DEBUGGING
    iSave = iSave + 1
    if (iSave > nSave) then
        exit
    end if !iSave > nSave
    do s=1, nSpecies
        currentTrajectory(s, iSave) = n(s)
        countSave = countSave + 1 !DEBUGGING
    end do !s=1, nSpecies
    if (iSave > 120) then !DEBUGGING
        if (currentTrajectory(1, iSave) == currentTrajectory(1, iSave-1) ) then
        end if
    end if !DEBUGGING
end do !while (iSave < int(t/dt))

end subroutine saveState






subroutine calcPropensities
use simulation
use model
implicit none

!This routine calculates the reaction propensities and their sum.
!It also finds the sum of critical reaction propensities, 
!and identifies the lists of critical and noncritical reactions.

real*8 pro !tmp value of a reaction propensity
integer :: rxn !interator for looping over reactions

criticalPropensitySum = 0.d0
propensitySum = 0.d0

do rxn = 1, nRxn
    pro = rates(rxn) * n( reactants(rxn) )
    propensities(rxn) = pro
    propensitySum = propensitySum + pro
    if( n( reactants(rxn) ) <= nCritical ) then
        isCritical(rxn) = .true.
        criticalPropensitySum = criticalPropensitySum + pro
    else
        isCritical(rxn) = .false.
    end if !n( reactants(rxn) ) <= nCritical
end do !rxn = 1, nRxn=

end subroutine calcPropensities






real*8 function tentativeLeap()
use simulation
use model
implicit none

real*8 :: pro, term, vrs !local tmp variables for efficiency
real*8 :: tmpTau, minTau, rns !real(n(s),kind=8) -- just a type converted population
integer :: i, rxn, s !iterators

!calculate max of N (this gives max of eps*N and of eps^2*N^2, for first order rxns
do s = 1, nSpecies
    maxEpsN(s) = eps * real(n(s),kind=8)
    if( maxEpsN(s) < 1.d0 ) then
        maxEpsN(s) = 1.d0
    end if ! maxEpsN(s) < 1.d0
end do !s = 1, nSpecies
mu(:) = 0.d0
sigmaSq(:) = 0.d0
do rxn = 1, nRxn
    if(.not.(isCritical(rxn))) then
        do s = 1, nSpecies
            vrs = v(rxn,s)
            pro = propensities(rxn)
            term = vrs * pro
            mu( s ) = mu( s ) + term
            sigmaSq( s ) = sigmaSq( s ) + vrs * term
        end do
    end if !critical(rxn)
end do !rxn = 1, nRxn
mu(:) = abs( mu(:) )
minTau = min( maxEpsN(1) / mu(1), maxEpsN(1) * maxEpsN(1)/ sigmaSq(1) )
if( nSpecies > 1 ) then
    do s = 2, nSpecies   
        tmpTau = min( maxEpsN(s) / mu(s), maxEpsN(s) * maxEpsN(s) / sigmaSq(s) )
        if( tmpTau < minTau ) then
            minTau = tmpTau
        end if !tmpTau < minTau
    end do !s = 1, nSpecies
end if !nSpecies > 1
tentativeLeap = minTau

end function tentativeLeap





!This function that samples from a Poisson distribution
!with parameter lambda = xm is adapted from
!Press, Teukolsky, Vetterling, and Flannery,
!Numerical Recipes in Fortran 77: The Art of Scientific Computing,
!Second Edition, Vol. 1, 1997, page 284
!This is implemented to return an integer. It can of course
!be modified to return the number as integer*8 if large populations
!are needed, or as a real, if working with the integers is inconvenient
!for some reason
function poidev(xm)
implicit none
real*8 xm,pi,rand
real*8 normalDev
integer*8 poidev
parameter (pi=3.14159265358979d0)
real*8 alxm,em,g,oldm,sq,t,y,check
save alxm,g,oldm,sq
data oldm /-1.d0/
if (xm .gt. 1.e6) then
    em = normalDev( xm, sqrt(xm) )
else if (xm .lt. 12.d0) then !If lambda is small, calculate by the "direct method"
  if (xm .ne. oldm) then !Check whether you happen to use same lambda twice
    oldm = xm      !If not, then recalculate exponential
    g = exp(-xm)
  end if !xm .ne. oldm
  em = -1
  t = 1.d0
  em = em + 1.d0
  call random_number(rand)
  t = t * rand
  do while (t .gt. g)
    em = em + 1.d0
    call random_number(rand)
    t = t * rand
  end do !while (t .gt. g)
else
!Otherwise, use rejection method
  if ( xm .ne. oldm ) then
    oldm = xm
    sq = sqrt(2.d0*xm)
    alxm = log(xm)
    g = xm * alxm - log_gamma( xm + 1.d0 )
  end if !xm .ne. oldm
  check = 2.d0
  t = 1.d0
!The code in numerical recipes uses goto statements in place of these
!while loops. The code could be made more compact, but probably there
!would be no significant change in performance.
  do while ( check .gt. t )
    call random_number(rand)
    y = tan( pi * rand )
    em = sq * y + xm
    do while ( em .lt. 0.d0 )
      call random_number(rand)
      y = tan( pi * rand )
      em = sq * y + xm
    end do !while ( em .lt. 0.d0 )
    em = int( em )
    t = 0.9d0 * (1.d0 + y * y) * exp( em * alxm - log_gamma( em + 1.d0 ) - g)
    call random_number(check)
  end do !while ( check .gt. t )
end if !(xm .lt. 12.d0)
poidev = nint(em,kind=8)
return
end function poidev

!normalDev returns a sample from a normal distribution.
!Uses the algorithm from Leva, JL 1992. "A Fast Normal Random Number Generator",
!ACM Trans. Math. Softw., 18 (4), pp 449-453.
real*8 function normalDev( mu, sig )
implicit none

  real*8 :: u, v, x, y , q, mu, sig
  logical :: reject
    reject = .true.
  do while( reject )
    call random_number(u)
    call random_number(v)
    v = 1.7156d0 * (v - 0.5)
    x = u - 0.449871d0
    y = abs(v) + 0.386595d0
    q = x * x + y * (0.19600d0 * y - 0.25472d0 * x)
    if (q .lt. 0.27597d0) then
      reject = .false.
    else if (q .le. 0.27846d0) then
      if (v * v .le. -4.d0 * log(u) * u * u) then
        reject = .false.
      end if !v * v < -4. * log(u) * u * u
    end if !(q .lt. 0.27846d0)
  end do !while( reject )
  normalDev = mu + sig * v / u

end function normalDev
