import pandas as pd

samples = pd.read_table(config['samples'], sep="\t").set_index("sample_id", drop=False)
units = pd.read_table(config['units'], sep="\t").set_index("sample_id", drop=False)

umi_bc_read = dict(zip(units['sample_id'], units['umi_bc_read']))
tx_read = dict(zip(units['sample_id'], units['tx_read']))
condition = dict(zip(samples['sample_id'], samples['condition']))