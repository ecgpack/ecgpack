#!/bin/python3
#
# Values of fundamental constants and energy conversion factors according to CODATA 2022
electron_mass=5.485799090441E-04    #electron mass in atomic mass units, u
d_electron_mass=00.000000000097E-04
hartree2u=2.92126231797E-08        #conversion factor hartree to atomic mass unit, u
d_hartree2u=0.00000000091E-08
ev2u=1.07354410083E-09             #conversion factor eV to atomic mass unit, u
d_ev2u=0.00000000033E-09

# Note: to get the nuclear mass accurate to 4 decimal figures, the total electron binding
# energy only needs to be accurate to about 0.1 hartree or even less.
# Ideally, the total electron binding energy should include the finite-nuclear-mass and
# relativistic effects. But, again, in reality one might ignore these small effects
# without consequences.

print('This script calculates masses of light nuclei (in atomic units, i.e. in the electron masses) including error estimates.') 
print('The values of the fundamental constants are taken from CODATA 2022 recommended values:') 
print('    P. Mohr, D. Newell, B. Taylor, and E. Tiesinga')
print('    CODATA Recommended Values of the Fundamental Physical Constants: 2022')
print('    arXiv:2409.03787 [hep-ph]')
print('    https://arxiv.org/abs/2409.03787')
print('Atomic weights are taken from AME 2020 atomic mass evaluation')
print('    Wang et al., Chin. Physics C 45, 030003 (2021)')
print('    https://dx.doi.org/10.1088/1674-1137/abddaf')
print('Most common isotopic compositions are avaialble at this NIST website')
print('    https://physics.nist.gov/cgi-bin/Compositions/stand_alone.pl')

print('')
print('Masses of most common nuclear isotopes:')

Z=1
atomic_mass=1007825.031898E-06
d_atomic_mass=0.000014E-06
electron_binding_energy=0.49973449614E0*hartree2u
d_electron_binding_energy=0.000002E0*hartree2u
nuclear_mass=(atomic_mass+electron_binding_energy)/electron_mass-Z
d_nuclear_mass=(d_atomic_mass+d_electron_binding_energy)/electron_mass+(atomic_mass+electron_binding_energy)*d_electron_mass/electron_mass**2
Mstr='M(1H)'
print('')
print('{:6s} = {:14.8f}'.format(Mstr,nuclear_mass))
print('{:12s}±{:10.8f}'.format('',d_nuclear_mass))
print('')
CODATAstr='1836.152673426(32)  [CODATA 2022]'
print('{:6s} =  '.format(Mstr) + CODATAstr)

Z=1
atomic_mass=2014101.777844E-06
d_atomic_mass=0.000015E-06
electron_binding_energy=0.49987047167E0*hartree2u
d_electron_binding_energy=0.000002E0*hartree2u
nuclear_mass=(atomic_mass+electron_binding_energy)/electron_mass-Z
d_nuclear_mass=(d_atomic_mass+d_electron_binding_energy)/electron_mass+(atomic_mass+electron_binding_energy)*d_electron_mass/electron_mass**2
Mstr='M(2H)'
print('')
print('{:6s} = {:14.8f}'.format(Mstr,nuclear_mass))
print('{:12s}±{:10.8f}'.format('',d_nuclear_mass))
print('')
CODATAstr='3670.482967655(63)  [CODATA 2022]'
print('{:6s} =  '.format(Mstr) + CODATAstr)

Z=1
atomic_mass=3016049.281320E-06
d_atomic_mass=0.000081E-06
electron_binding_energy=0.49991571296E0*hartree2u
d_electron_binding_energy=0.000002E0*hartree2u
nuclear_mass=(atomic_mass+electron_binding_energy)/electron_mass-Z
d_nuclear_mass=(d_atomic_mass+d_electron_binding_energy)/electron_mass+(atomic_mass+electron_binding_energy)*d_electron_mass/electron_mass**2
Mstr='M(3H)'
print('')
print('{:6s} = {:14.8f}'.format(Mstr,nuclear_mass))
print('{:12s}±{:10.8f}'.format('',d_nuclear_mass))
print('')
CODATAstr='5496.92153551(21)  [CODATA 2022]'
print('{:6s} =  '.format(Mstr) + CODATAstr)

Z=2
atomic_mass=3016029.321967E-06
d_atomic_mass=0.000060E-06
electron_binding_energy=2.903271162695E0*hartree2u #includes rel. corr., source: M. Stanke et al. JCP 126, 194312 (2007)
d_electron_binding_energy=0.00002E0*hartree2u
nuclear_mass=(atomic_mass+electron_binding_energy)/electron_mass-Z
d_nuclear_mass=(d_atomic_mass+d_electron_binding_energy)/electron_mass+(atomic_mass+electron_binding_energy)*d_electron_mass/electron_mass**2
Mstr='M(3He)'
print('')
print('{:6s} = {:14.8f}'.format(Mstr,nuclear_mass))
print('{:12s}±{:10.8f}'.format('',d_nuclear_mass))
print('')
CODATAstr='5495.88527984(16)  [CODATA 2022]'
print('{:6s} =  '.format(Mstr) + CODATAstr)

Z=2
atomic_mass=4002603.254130E-06
d_atomic_mass=0.000158E-06
electron_binding_energy=2.903408504575E0*hartree2u #includes rel. corr., source: M. Stanke et al. JCP 126, 194312 (2007)
d_electron_binding_energy=0.00002E0*hartree2u
nuclear_mass=(atomic_mass+electron_binding_energy)/electron_mass-Z
d_nuclear_mass=(d_atomic_mass+d_electron_binding_energy)/electron_mass+(atomic_mass+electron_binding_energy)*d_electron_mass/electron_mass**2
Mstr='M(4He)'
print('')
print('{:6s} = {:14.8f}'.format(Mstr,nuclear_mass))
print('{:12s}±{:10.8f}'.format('',d_nuclear_mass))
print('')
CODATAstr='7294.29954171(17)  [CODATA 2022]'
print('{:6s} =  '.format(Mstr) + CODATAstr)

Z=3
atomic_mass=6015122.8874E-06
d_atomic_mass=0.0015E-06
electron_binding_energy=7.47799240823E0*hartree2u #includes rel. corr., source:  M. Stanke et. al., PRA 78, 052507 (2008)
d_electron_binding_energy=0.0001E0*hartree2u
nuclear_mass=(atomic_mass+electron_binding_energy)/electron_mass-Z
d_nuclear_mass=(d_atomic_mass+d_electron_binding_energy)/electron_mass+(atomic_mass+electron_binding_energy)*d_electron_mass/electron_mass**2
Mstr='M(6Li)'
print('')
print('{:6s} = {:14.8f}'.format(Mstr,nuclear_mass))
print('{:12s}±{:10.8f}'.format('',d_nuclear_mass))

Z=3
atomic_mass=7016003.434E-06
d_atomic_mass=0.004E-06
electron_binding_energy=7.47809366029E0*hartree2u #includes rel. corr., source:  M. Stanke et. al., PRA 78, 052507 (2008)
d_electron_binding_energy=0.0001E0*hartree2u
nuclear_mass=(atomic_mass+electron_binding_energy)/electron_mass-Z
d_nuclear_mass=(d_atomic_mass+d_electron_binding_energy)/electron_mass+(atomic_mass+electron_binding_energy)*d_electron_mass/electron_mass**2
Mstr='M(7Li)'
print('')
print('{:6s} = {:14.8f}'.format(Mstr,nuclear_mass))
print('{:12s}±{:10.8f}'.format('',d_nuclear_mass))

Z=4
atomic_mass=9012183.06E-06
d_atomic_mass=0.08E-06
electron_binding_energy=14.668795768E0*hartree2u #includes rel. corr., source:  M. Stanke et. al., PRA 80, 022514 (2009)
d_electron_binding_energy= 0.0005E0*hartree2u
nuclear_mass=(atomic_mass+electron_binding_energy)/electron_mass-Z
d_nuclear_mass=(d_atomic_mass+d_electron_binding_energy)/electron_mass+(atomic_mass+electron_binding_energy)*d_electron_mass/electron_mass**2
Mstr='M(9Be)'
print('')
print('{:6s} = {:14.8f}'.format(Mstr,nuclear_mass))
print('{:12s}±{:10.8f}'.format('',d_nuclear_mass))

Z=5
atomic_mass=10012936.862E-06
d_atomic_mass=0.016E-06
electron_binding_energy=(24.65250024E0+0.005E0)*hartree2u #nonrel. energy:  S. Bubin and L. Adamowicz, PRA 83, 022505 (2011); rel. corr. P. Seth et. al., JCP 134, 084105 (2011)
d_electron_binding_energy= 0.001E0*hartree2u
nuclear_mass=(atomic_mass+electron_binding_energy)/electron_mass-Z
d_nuclear_mass=(d_atomic_mass+d_electron_binding_energy)/electron_mass+(atomic_mass+electron_binding_energy)*d_electron_mass/electron_mass**2
Mstr='M(10B)'
print('')
print('{:6s} = {:14.8f}'.format(Mstr,nuclear_mass))
print('{:12s}±{:10.8f}'.format('',d_nuclear_mass))

Z=5
atomic_mass=11009305.167E-06
d_atomic_mass=0.013E-06
electron_binding_energy=(24.65262387E0+0.005E0)*hartree2u #nonrel. energy:  S. Bubin and L. Adamowicz, PRA 83, 022505 (2011); rel. corr. P. Seth et. al., JCP 134, 084105 (2011)
d_electron_binding_energy= 0.001E0*hartree2u
nuclear_mass=(atomic_mass+electron_binding_energy)/electron_mass-Z
d_nuclear_mass=(d_atomic_mass+d_electron_binding_energy)/electron_mass+(atomic_mass+electron_binding_energy)*d_electron_mass/electron_mass**2
Mstr='M(11B)'
print('')
print('{:6s} = {:14.8f}'.format(Mstr,nuclear_mass))
print('{:12s}±{:10.8f}'.format('',d_nuclear_mass))

Z=6
atomic_mass=12.0E0
d_atomic_mass=0.000000000E0
electron_binding_energy=(37.84446E0-0.0017166E0+0.015E0)*hartree2u
d_electron_binding_energy= 0.005E0*hartree2u
nuclear_mass=(atomic_mass+electron_binding_energy)/electron_mass-Z
d_nuclear_mass=(d_atomic_mass+d_electron_binding_energy)/electron_mass+(atomic_mass+electron_binding_energy)*d_electron_mass/electron_mass**2
Mstr='M(12C)'
print('')
print('{:6s} = {:14.8f}'.format(Mstr,nuclear_mass))
print('{:12s}±{:10.8f}'.format('',d_nuclear_mass))

Z=6
atomic_mass=13003354.83534E-06
d_atomic_mass=0.00025E-06
electron_binding_energy=(37.84446E0-0.00158412E0+0.015E0)*hartree2u
d_electron_binding_energy= 0.005E0*hartree2u
nuclear_mass=(atomic_mass+electron_binding_energy)/electron_mass-Z
d_nuclear_mass=(d_atomic_mass+d_electron_binding_energy)/electron_mass+(atomic_mass+electron_binding_energy)*d_electron_mass/electron_mass**2
Mstr='M(13C)'
print('')
print('{:6s} = {:14.8f}'.format(Mstr,nuclear_mass))
print('{:12s}±{:10.8f}'.format('',d_nuclear_mass))

Z=6
atomic_mass=14003241.989E-06
d_atomic_mass=0.004E-06
electron_binding_energy=(37.84446E0-0.00147098E0+0.015E0)*hartree2u
d_electron_binding_energy=0.005E0*hartree2u
nuclear_mass=(atomic_mass+electron_binding_energy)/electron_mass-Z
d_nuclear_mass=(d_atomic_mass+d_electron_binding_energy)/electron_mass+(atomic_mass+electron_binding_energy)*d_electron_mass/electron_mass**2
Mstr='M(14C)'
print('')
print('{:6s} = {:14.8f}'.format(Mstr,nuclear_mass))
print('{:12s}±{:10.8f}'.format('',d_nuclear_mass))

Z=7
atomic_mass=14003074.00425E-06
d_atomic_mass=0.00024E-06
electron_binding_energy=(54.58867E0+0.03E0)*hartree2u
d_electron_binding_energy= 0.01E0*hartree2u
nuclear_mass=(atomic_mass+electron_binding_energy)/electron_mass-Z
d_nuclear_mass=(d_atomic_mass+d_electron_binding_energy)/electron_mass+(atomic_mass+electron_binding_energy)*d_electron_mass/electron_mass**2
Mstr='M(14N)'
print('')
print('{:6s} = {:14.8f}'.format(Mstr,nuclear_mass))
print('{:12s}±{:10.8f}'.format('',d_nuclear_mass))

Z=7
atomic_mass=15000108.8983E-06
d_atomic_mass=0.0006E-06
electron_binding_energy=(54.58867E0+0.03E0)*hartree2u
d_electron_binding_energy= 0.01E0*hartree2u
nuclear_mass=(atomic_mass+electron_binding_energy)/electron_mass-Z
d_nuclear_mass=(d_atomic_mass+d_electron_binding_energy)/electron_mass+(atomic_mass+electron_binding_energy)*d_electron_mass/electron_mass**2
Mstr='M(15N)'
print('')
print('{:6s} = {:14.8f}'.format(Mstr,nuclear_mass))
print('{:12s}±{:10.8f}'.format('',d_nuclear_mass))

Z=8
atomic_mass=15994914.6193E-06
d_atomic_mass=0.0003E-06
electron_binding_energy=(75.0673E0+0.053E0)*hartree2u
d_electron_binding_energy= 0.02E0*hartree2u
nuclear_mass=(atomic_mass+electron_binding_energy)/electron_mass-Z
d_nuclear_mass=(d_atomic_mass+d_electron_binding_energy)/electron_mass+(atomic_mass+electron_binding_energy)*d_electron_mass/electron_mass**2
Mstr='M(16O)'
print('')
print('{:6s} = {:14.8f}'.format(Mstr,nuclear_mass))
print('{:12s}±{:10.8f}'.format('',d_nuclear_mass))

Z=8
atomic_mass=16999131.7560E-06
d_atomic_mass=0.0007E-06
electron_binding_energy=(75.0673E0+0.053E0)*hartree2u
d_electron_binding_energy= 0.02E0*hartree2u
nuclear_mass=(atomic_mass+electron_binding_energy)/electron_mass-Z
d_nuclear_mass=(d_atomic_mass+d_electron_binding_energy)/electron_mass+(atomic_mass+electron_binding_energy)*d_electron_mass/electron_mass**2
Mstr='M(17O)'
print('')
print('{:6s} = {:14.8f}'.format(Mstr,nuclear_mass))
print('{:12s}±{:10.8f}'.format('',d_nuclear_mass))

Z=8
atomic_mass=17999159.6121E-06
d_atomic_mass=0.0007E-06
electron_binding_energy=(75.0673E0+0.053E0)*hartree2u
d_electron_binding_energy= 0.02E0*hartree2u
nuclear_mass=(atomic_mass+electron_binding_energy)/electron_mass-Z
d_nuclear_mass=(d_atomic_mass+d_electron_binding_energy)/electron_mass+(atomic_mass+electron_binding_energy)*d_electron_mass/electron_mass**2
Mstr='M(18O)'
print('')
print('{:6s} = {:14.8f}'.format(Mstr,nuclear_mass))
print('{:12s}±{:10.8f}'.format('',d_nuclear_mass))

print('')
print('Masses of some exotic particles:')

print('')
Mstr='M(muon)'
CODATAstr='206.7682827(46)  [CODATA 2022]'
print('{:8s} = '.format(Mstr) + CODATAstr)

print('')
Mstr='M(tauon)'
CODATAstr='3477.23(23)  [CODATA 2022]'
print('{:8s} = '.format(Mstr) + CODATAstr)

print('')
Mstr='M(pion±)'
CODATAstr='273.13287(47)  [Phys. Lett. B 796, 11 (2019)]'  #139.57061±0.00024 MeV/c^2 From M. Daum, R. Frosch 1, P.-R. Kettle, Phys. Lett. B 796, 11 (2019)
print('{:8s} = '.format(Mstr) + CODATAstr)
