#!/usr/bin/env python3

import argparse
import datetime
import sys

import polars as pl

parser = argparse.ArgumentParser()
parser.add_argument('outprefix')
parser.add_argument('samples')
parser.add_argument('data_fname')
parser.add_argument('phenotype_name')
parser.add_argument('fam_file')
parser.add_argument('pcs_fname')
parser.add_argument('n_pcs', type=int)
parser.add_argument('assessment_ages')
parser.add_argument(
    '--categorical-covar-names',
    nargs='*',
    default=[]
)
parser.add_argument(
    '--categorical-covar-files',
    nargs='*',
    default=[]
)
args = parser.parse_args()

assert 0 <= args.n_pcs and args.n_pcs <= 40
assert len(args.categorical_covar_names) == len(args.categorical_covar_files)

cat_drop_num = 50
cat_drop_frac = 0.001 #0.1%

df = pl.read_csv(
    args.samples,
    new_columns = ['id']
).select([
    pl.col('id')
]).join(
    pl.read_csv(
        args.data_fname,
        separator = '\t',
        new_columns = ['id', args.phenotype_name]
    ).select([
       pl.col('id'), pl.col(args.phenotype_name)
    ]),
    how = 'inner',
    on = 'id'
).join(
    pl.read_csv(
        args.fam_file,
        separator = ' ',
        has_header = False,
        new_columns = ['id', 'id_duplicate', 'zero', 'zero_duplicate', 'sex', 'batch'],
    ).select([
        pl.col('id'), pl.col('sex')
    ]),
    how = 'inner',
    on = 'id'
).join(
    pl.read_csv(
        args.pcs_fname,
        separator = '\t',
        new_columns = ['id', *[f'pc_{n}' for n in range(1, 41)]]
    ).select([
        pl.col(col) for col in
        ['id', *[f'pc_{n}' for n in range(1, args.n_pcs+1)]]
    ]),
    how = 'inner',
    on = 'id'
).join(
    pl.read_csv(
        args.assessment_ages,
        separator = '\t',
        new_columns = ['id', 'init-assessment-age', 'repeat-assessment-age', 'imaging-visit-age', 'repeat-imaging-age'],
    ).select([
        pl.col('id'), pl.col('init-assessment-age').alias('age')
    ]).filter(
        ~pl.col('age').is_null()
    ),
    how = 'inner',
    on = 'id'
)

for (f, covar_name) in zip(args.categorical_covar_files, args.categorical_covar_names):
    df = df.join(
        pl.read_csv(
            f,
            separator='\t',
            new_columns = ['id', covar_name] # assuming second column is the value for this covar for the first assessment, following columns will be skipped
        ).select([
            pl.col('id'), pl.col(covar_name)    
        ]),
        how = 'inner',
        on = 'id'
    ).filter(
        ~pl.col(covar_name).is_null()
    )

n_samples = df.shape[0]
filter_condition = None
for covar in args.categorical_covar_names:
    grouping = df.group_by(covar).agg(pl.count())
    vals_to_counts = dict(zip(grouping[covar], grouping['count']))

    # remove the most common covariate value
    for val in vals_to_counts:
        if vals_to_counts[val] == max(vals_to_counts.values()):
            break
    print(f'Breakdown for categorical covar {covar}: {vals_to_counts}')
    print(f'Using {val} as the default value')
    del vals_to_counts[val]

    for val in vals_to_counts:
        n_with_val = df.filter(pl.col(covar) == val).shape[0]
        if n_with_val <= cat_drop_num or n_with_val <= cat_drop_frac*n_samples:
            print(f'Dropping samples with covar {covar} value {val}')
            new_condition = pl.col(covar) != val
            if filter_condition is None:
                filter_condition = new_condition
            else:
                filter_condition = filter_condition & new_condition
        else:
            df = df.with_columns((pl.col(covar) == val).cast(int).alias(f'{covar}_is_{val}'))
    df = df.drop(covar)

if filter_condition is not None:
    df = df.filter(filter_condition)

df.write_csv(f'{args.outprefix}.tab', separator='\t')

