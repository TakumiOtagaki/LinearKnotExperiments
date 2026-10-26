#!/usr/bin/env python3
"""Relocatable adapter for the fixed resource comparison (measured wrapper)."""
import argparse
import os
from pathlib import Path


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('-i', type=Path, required=True)
    parser.add_argument('--method', required=True,
                        choices=('cparty', 'exact', 'beam100', 'spark', 'fold100'))
    args = parser.parse_args()
    build = Path(__file__).resolve().parent
    parameter = os.environ['CPARTY_FIXED_ENERGY_PARAM_FILE']
    if args.method == 'cparty':
        command = [str(build/'E0CPartyOracle'), '-i', str(args.i)]
    elif args.method in ('exact', 'beam100'):
        command = [str(build/'PFScaling'), '-i', str(args.i), '--method', args.method]
    elif args.method == 'fold100':
        command = [str(build/'linearknot'), 'fold', '-i', str(args.i), '-b', '100',
                   '--max-hairpin', '100', '--max-case3-depth', '100', '--ranking', 'outer']
    else:
        _, sequence, scaffold = args.i.read_text().splitlines()[:3]
        command = [str(build/'Spark'), '-d', '2', '-P', parameter, '-r', scaffold, sequence]
    os.execv(command[0], command)


if __name__ == '__main__':
    main()
