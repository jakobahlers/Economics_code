import pandas as pd
import numpy as np
from tqdm import tqdm
import re
import os
from helpers.sql import read_db, connect_to_db, vacuum_db, update_database

engine = connect_to_db()

snp = read_db(statement='select distinct(ric) from snp', engine=engine, idx_col=None)
snp.head()

PATH = 'XXXXX'

def format_dataframe(df):
    # indexing by date
    df = df.set_index('Date', drop=True)
    df.index.name = None

    # there can be duplicated timestamps - needs to be removed
    df = df[~df.index.duplicated(keep='last')]

    # set correct names
    df = df.rename(columns={'Total Return': df['Instrument'].unique()[0]})

    # drop Instrument column
    df = df.drop(columns=['Instrument'])

    return df


# create placeholder DataFrame
daily_returns = pd.DataFrame(index=pd.date_range(
    start='1995-01-01', end='2023-02-01'))


for ric in tqdm(snp['ric']):
    eikon = pd.read_excel(f'{PATH}{ric}_daily.xlsx', index_col=0)

    try:
        eikon = format_dataframe(df=eikon)

        # merge data
        daily_returns = daily_returns.merge(
            eikon, how='left', left_index=True, right_index=True)
    except:
        print(ric)


# The following tickers have 0 observations. 
# However they've exited the index before the analysis starts
for col in daily_returns.columns:
    if col not in daily_returns.dropna(axis=1, how='all').columns:
        print(f'\'{col}\',')


# drop from the data-set
daily_returns = daily_returns.dropna(axis=1,how='all')


# drop rows with no observations (i.e. weekends, holidays...)
daily_returns = daily_returns.dropna(axis=0, how='all')

# div by 100
daily_returns = daily_returns / 100


# use regular incremental index
daily_returns['date'] = daily_returns.index
daily_returns = daily_returns.reset_index(drop=True)

# move date column to beginning of 
first_column = daily_returns.pop('date')
daily_returns.insert(0, 'date', first_column)

daily_returns.to_sql(name='returns', con=engine, if_exists='replace')
vacuum_db(engine=engine)
update_database()

