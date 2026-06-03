#!/usr/bin/env python3

## Read outputs per workflow rule table and print a table in markdown format for the specified
## analysis. Used to generate the analysis product tables shown in the README.md file

import sys
import csv
import argparse

def tsv_to_markdown(filepath, analysis=None):
    with open(filepath, newline='') as f:
        reader = csv.DictReader(f, delimiter='\t')
        rows = list(reader)

    if not rows:
        print("Empty file.")
        return

    if analysis:
        rows = [r for r in rows if r.get('Analysis') == analysis]
        if not rows:
            print(f"No rows found for analysis: {analysis}")
            return

    output_cols = [col for col in rows[0].keys() if col != 'Analysis']

    print('| ' + ' | '.join(output_cols) + ' |')
    print('| ' + ' | '.join(['---'] * len(output_cols)) + ' |')
    for row in rows:
        cells = []
        for i, col in enumerate(output_cols):
            cells.append(f'`{row[col]}`' if i == 1 else row[col])
        print('| ' + ' | '.join(cells) + ' |')

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Convert a TSV file to a GitHub markdown table.')
    parser.add_argument('file', help='Path to the .tsv file')
    parser.add_argument('--analysis', help='Filter rows by the Analysis column')
    args = parser.parse_args()

    tsv_to_markdown(args.file, args.analysis)
    
