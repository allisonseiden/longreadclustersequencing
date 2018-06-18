import pandas as pd
import numpy as np

dates = pd.date_range('20130101', periods=6)
#print(dates)

df = pd.DataFrame(np.random.randn(6,4), index=dates, columns=list('ABCD'))
#print(df)

test = pd.read_table('~/Downloads/1-00801.hg38.dnv.bed')
print(test)
