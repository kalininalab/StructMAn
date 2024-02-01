#!/usr/bin/python

# Copyright 2010 Max-Planck-Institut Informatik
#
#    This file is part of RINerator
#
#    RINerator is free software: you can redistribute it and/or modify
#    it under the terms of the GNU General Public License as published by
#    the Free Software Foundation, either version 3 of the License, or
#    (at your option) any later version.
#
#    RINerator is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#    GNU General Public License for more details.
#
#    You should have received a copy of the GNU General Public License
#    along with RINerator.  If not, see <http://www.gnu.org/licenses/>.

import sys

from . import nc_interactions, nci_graph, nci_residues, prt_set, st_selection

if len(sys.argv) != 2:
    print("USAGE: get_ncint.py script")

s_file = sys.argv[1]
print("INFO: Executing %s" % s_file)

execfile(s_file)
