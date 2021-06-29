# this is a python script which shows how to read the waveburst root file

# load ROOT stuff
from ROOT import *

# load HEALPix stuff
gSystem.Load("$HOME_CFITSIO/libcfitsio.so")
gSystem.Load("$HOME_HEALPIX/src/cxx/cxxsupport/lib/libcxxsupport.so")
gSystem.Load("$HOME_HEALPIX/src/cxx/libfftpack/lib/libfftpack.so")
gSystem.Load("$HOME_HEALPIX/src/cxx/c_utils/lib/libc_utils.so")
gSystem.Load("$HOME_HEALPIX/src/cxx/libpsht/lib/libpsht.so")
gSystem.Load("$HOME_HEALPIX/src/cxx/Healpix_cxx/lib/HEALPix.so")

# load wavelet library
gSystem.Load("wavelet.so")

# open waveburst root file
f = TFile( 'wave_file.root' )
# get waveburst tree 
tree = f.Get("waveburst")
# define tree cuts
treeformula = TTreeFormula("cuts", "rho[1]>5 && netcc[0]>0.8", tree);

# define histogram for netcc[0] plot
hist = TH1F( 'hist', 'netcc', 100, 0., 1. )

# get the number of entries
size = tree.GetEntries()
print size

# loop over the tree entries
for j in range(0,size-1):
  # get entry j
  tree.GetEntry(j)
  # select entry
  if treeformula.EvalInstance() != 0: 
    print tree.netcc[0]
    print tree.rho[1]
    # fill histogram
    hist.Fill(tree.netcc[0])

# plot and save histogram
c = TCanvas()
hist.Draw()
c.Print("netcc0.png")

# plot and save rho[1]
tree.Draw("rho[1]","rho[1]>5 && netcc[0]>0.8")
c.Print("rho1.png")

# plot and save rho[1] netcc[0]
tree.Draw("rho[1]:netcc[0]","rho[1]>5 && netcc[0]>0.8","colz")
c.Print("rho1_vs_netcc0.png")

