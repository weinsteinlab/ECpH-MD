 package require psfgen 	 
 topology top_all36_prot.rtf
 topology toppar_water_ions.str
 segment HEW { pdb 3lzt-prot.pdb; first NTER; last CTER } 	 
 coordpdb 3lzt-prot.pdb HEW
 patch GLUP HEW:7
 patch GLUP HEW:35
 patch ASPP HEW:18
 patch ASPP HEW:48
 patch ASPP HEW:52
 patch ASPP HEW:66
 patch ASPP HEW:87
 patch ASPP HEW:101
 patch ASPP HEW:119
 patch DISU HEW:30 HEW:115
 patch DISU HEW:6 HEW:127
 patch DISU HEW:64 HEW:80
 patch DISU HEW:76 HEW:94
 segment CRW { pdb  3lzt-hoh.pdb }
 coordpdb 3lzt-hoh.pdb CRW
# segment ION { pdb  ions.pdb }
# coordpdb ions.pdb ION
 regenerate angles dihedrals
 guesscoord
 writepdb hewl-crw.pdb
 writepsf hewl-crw.psf
