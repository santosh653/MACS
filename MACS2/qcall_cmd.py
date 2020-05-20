# Time-stamp: <2019-11-04 10:41:13 taoliu>

"""Description: MACS 2 main executable

This code is free software; you can redistribute it and/or modify it
under the terms of the BSD License (see the file LICENSE included with
the distribution).
"""

# ------------------------------------
# python modules
# ------------------------------------

import os
import sys
import logging
from time import strftime
import tempfile

# ------------------------------------
# own python modules
# ------------------------------------
from MACS2.OptValidator import opt_validate_qcall as opt_validate
from MACS2.Pileup import unified_pileup_bdg
from MACS2.Constants import *
# ------------------------------------
# Main function
# ------------------------------------
def run( args ):
    """The Quick Peak Calling
    
    """
    # Parse options...
    options = opt_validate( args )
    # end of parsing commandline options
    info = options.info
    warn = options.warn
    debug = options.debug
    error = options.error
    #0 output arguments
    options.PE_MODE = options.format in ('BAMPE','BEDPE')

    #tempfile.tempdir = options.tempdir

    #1 Read tag files
    info("# read alignment files...")
    if options.PE_MODE:
        info("# read input file in Paired-end mode.")
        (treat, extsize) = load_frag_files_options ( options ) # return PETrackI object
        t0 = treat.total # total fragments
        info("# total fragments/pairs in alignment file: %d" % (t0) )
        info("# Pileup paired-end alignment file.")
        btrack = unified_pileup_bdg(treat, None , 1, halfextension=False)        
    else:
        (treat, tsize) = load_tag_files_options  (options)
        t0 = treat.total
        info("# total tags in alignment file: %d", t0)
        info("# Pileup alignment file, extend each read towards downstream direction with %d bps" % options.extsize)
        extsize = options.extsize
        btrack = unified_pileup_bdg(treat, options.extsize, 1, halfextension=False)

    info("Effective Genome Size(bp): %d" % options.gsize)
    genome_bg = extsize * t0 / options.gsize
    info("Average genome coverage: %.4f" % genome_bg)    

    coverage_cutoff = options.cutoff * genome_bg
    info("Call peaks from bedGraph with fold change cutoff %.2f, or the coverage cutoff %.2f" % (options.cutoff, coverage_cutoff))
    peaks = btrack.call_peaks(cutoff=float(coverage_cutoff),min_length=int(options.minlen),max_gap=int(options.maxgap),call_summits=False)

    info("Write peaks...")
    if options.ofile:
        options.oprefix = options.ofile
        nf = open( os.path.join( options.outdir, options.ofile ), 'w' )
    else:
        nf = open ( os.path.join( options.outdir, "%s_c%.1f_l%d_g%d_peaks.narrowPeak" % (options.oprefix,options.cutoff,options.minlen,options.maxgap)), "w" )
    peaks.write_to_narrowPeak(nf, name=options.oprefix.encode(), name_prefix=(options.oprefix+"_narrowPeak").encode(), score_column="score", trackline=False)
    info("Done")


def load_frag_files_options ( options ):
    """From the options, load treatment fragments and control fragments (if available).

    """
    options.info("#1 read treatment fragments...")

    tp = options.parser(options.ifile[0], buffer_size=options.buffer_size)
    treat = tp.build_petrack()
    tsize = tp.d

    if len(options.ifile) > 1:
        # multiple input
        for ifile in options.ifile[1:]:
            tp = options.parser(ifile, buffer_size=options.buffer_size)
            treat = tp.append_petrack( treat )

    treat.finalize()
    return (treat, tsize)

def load_tag_files_options ( options ):
    """From the options, load treatment tags and control tags (if available).

    """
    options.info("#1 read treatment tags...")
    tp = options.parser(options.ifile[0], buffer_size=options.buffer_size)
    tsize = tp.tsize()
    treat = tp.build_fwtrack()
    
    if len(options.ifile) > 1:
        # multiple input
        for ifile in options.ifile[1:]:
            tp = options.parser(ifile, buffer_size=options.buffer_size)
            treat = tp.append_fwtrack( treat )

    treat.finalize()
    return (treat, tsize)
