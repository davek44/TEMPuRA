#!/usr/bin/env python
from optparse import OptionParser
import os

################################################################################
# tempura.py
#
# Shared methods and data.
################################################################################

src_dir = os.environ['TEMPURADIR']

# data directories
dfam_dir = os.environ['DFAM']

# R script directory
r_dir = '%s/r' % src_dir
