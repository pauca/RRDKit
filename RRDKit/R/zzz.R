# This file is part of RRDKitInchi.
# 
# RRDKit is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
# 
# Foobar is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
# 
# You should have received a copy of the GNU General Public License
# along with RRDKitInchi.  If not, see <http://www.gnu.org/licenses/>.

#loadModule("p_smiles2mol",T) 
#loadModule("RWMol",T) 


.onAttach <- function(libname, pkgname) {
  ver <- read.dcf(file=system.file("DESCRIPTION", package=pkgname),
                  fields="Version")
  packageStartupMessage(paste(pkgname, ver))
  packageStartupMessage(paste("Using RDKit from $RDBASE =", Sys.getenv("RDBASE"),collapse = " "))
}
