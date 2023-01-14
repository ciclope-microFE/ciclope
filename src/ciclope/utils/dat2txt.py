#!/usr/bin/python
"""
This script extracts data from ccx dat-files.
Extended from the dat2txt.py helper script by
https://github.com/mkraska/CalculiX-Examples/tree/master/Scripts

"""
import re
import glob
import argparse
import logging
import textwrap

def dat2txt(filename, verbose=False):
    """Read CalculiX .DAT output file and write analysis results to .TXT file.

    Parameters
    ----------
    filename
        CalculiX ouput .DAT file.
    verbose
        Verbose output.

    Returns
    -------
    resname
        results name
    group
        results group name
    time_points
        number of time steps
    """

    if verbose:
        logging.basicConfig(level=logging.INFO)

    # variables initialization
    resname = ''
    group = ''
    data = {}
    pH = re.compile(' (.+) for .*set\\s(\\S+) and time  (.+)')
    skip = 0  # if empty lines are to be skipped
    body = 0  # if data lines are expected
    nev = 0  # number of eigenvalue file
    stop_on_empty = 0  # empty line triggers end of body
    time_points = 0

    # open job file
    logging.info('Reading file {}.'.format(filename))
    f = open(filename, "r")

    for line in f:
        if skip:  # if the line is known to be useless
            skip = skip - 1
            continue
        line = line.replace("(", "")
        line = line.replace(")", "")
        line = line.replace("\n", "")
        b = pH.match(line)
        if b:  # a result header was found
            resname = b.group(1)
            group = b.group(2)
            # print(group)
            res = resname + "_" + group
            time = float(b.group(3))
            if not (res in data.keys()):  # new result type
                res_filename = res + ".txt"
                data[res] = open(res_filename, "w")
            data[res].write("\n" + str(time) + " ")
            line_end = " "
            body = 1
            stop_on_empty = 0
        elif line.startswith("     E I G E N V A L U E   O U T P U T"):
            # eigenvalue data
            nev += 1
            res = "Eigenvalues_" + str(nev)
            # print(res)
            data[res] = open(res + ".txt", "w")
            data[res].write("# mode  EV  Re(f)_rad/time Re(f)_cycles/time Im(f)_rad/time\n")
            line_end = "\n"
            body = 1
            skip = 5
            stop_on_empty = 1
        elif line.startswith("     P A R T I C I P A T I O N   F A C T O R S"):
            # eigenvalue data
            res = "Eigenvalues_PF_" + str(nev)
            print(res)
            data[res] = open(res + ".txt", "w")
            data[res].write("# mode  UX  UY UZ RX RY RZ\n")
            line_end = "\n"
            body = 1
            skip = 3
            stop_on_empty = 1
        elif line.startswith("     E F F E C T I V E   M O D A L   M A S S"):
            # eigenvalue data
            res = "Eigenvalues_MM_" + str(nev)
            print(res)
            data[res] = open(res + ".txt", "w")
            data[res].write("# mode  UX  UY UZ RX RY RZ\n")
            line_end = "\n"
            body = 1
            skip = 3
            stop_on_empty = 1
        else:
            if body == 1:
                if not "force" in line and not "center" in line:
                    data[res].write(line + line_end)
                    if line.strip():
                        time_points += 1
                if stop_on_empty and line == "":
                    body = 0

    for name in data.keys():
        data[name].close()

    logging.info('Written {} for group {} to file: {}.\nTotal time points: {}.\n'.format(resname, group, res_filename, time_points))
    return resname, group, time_points

def main():
    description = textwrap.dedent('''\
                Script to process CalculiX Finite Element (FE) output .DAT files.
                Reads CalculiX FE analysis results and writes a .TXT file of the output variables and time steps.
                ''')

    parser = argparse.ArgumentParser(description=description, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('job', type=str, help='Input CalculiX job name.')
    # parser.add_argument('j', '--jobname', type=str, default=None, help='<CalculiX job name.')
    # parser.add_argument('fileout', type=str, default=None, help='<Required> Output filename (Abaqus .INP).')
    parser.add_argument('-v', '--verbose', dest='verbose', action='store_true', help='Verbose output.')
    parser.set_defaults(job=None, verbose=True)

    args = parser.parse_args()

    if args.verbose:
        logging.basicConfig(level=logging.INFO)

    # processing command line arguments, get the jobname
    if args.job is None:
        logging.info("No jobname given.\n")
        files = glob.glob("*.dat")
        if len(files) == 1:
            logging.info('File {} found'.format(files[0]))
            job = files[0]
        else:
            logging.exception('Available data files: {}'.format(*f) for f in files)
            quit()
    else:
        job = args.job + ".dat"

    # run dat2txt
    resname, group, time_points = dat2txt(job)

if __name__ == '__main__':
    main()

# names = ['time', 'fx', 'fy', 'fz']
# filename = '/home/gianthk/PycharmProjects/CT2FE/test_data/steel_foam/total force fx,fy,fz_NODES_B.txt'
# data2 = pd.read_table(filename, delim_whitespace=True, names=names, index_col=0)
