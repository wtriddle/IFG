""" Top level script run with py RunIFG.py [args, [...]] """

import sys
import getopt
# from IFG.src import smiles


def psuedoFunc(a, p):
    print("Allsheet = " + str(a))
    print("preciseSheet = " + str(p))


def main(argv):
    preciseSheet = False
    allSheet = False
    verbose = False
    try:
        args, opts = getopt.getopt(
            argv, 'aphv', ["all", "precise", "help", "verbose"])
    except getopt.GetoptError:
        print('Bad usage. Type smiles.py -h for help')
        sys.exit()
    for opt, arg in args:
        if opt in ['-h', '--help']:
            # Help option
            print(
                """smiles.py [ -a | --all <bool> ] [ -p | --precise <bool> ] [ -v || --verbose <bool> ] [-h | --help]
                """)
            sys.exit()
        elif opt in ['-a', '--all']:
            allSheet = True
        elif opt in ['-p', '--precise']:
            preciseSheet = True
        elif opt in ['-v', '--verbose']:
            verbose = True
    if not args:
        preciseSheet = True
        allSheet = True
    psuedoFunc(allSheet, preciseSheet)


if __name__ == '__main__':
    main(sys.argv[1:])
