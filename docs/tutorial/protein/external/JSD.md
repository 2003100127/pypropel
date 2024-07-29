## 1. Download the JSD package

The standalone package of JSD can be downloaded at [https://compbio.cs.princeton.edu/conservation/](https://compbio.cs.princeton.edu/conservation/). The JSD package uses the Jensen-Shannon divergence for computing residue conservation based on a MSA. 


!!! tip

    The JSD output can be used as a feature to describe residue evolutionary profile at each MSA column.


## 2. Running the JSD package

We wrote a wrapper to use the package in computer clusters for batch processing of the MSA files. It takes as input MSAs in `clustal` format.

Before using, we need to set a few parameters.

``` py linenums="1"
param_config = {
    'method': '-s',
    # 'window': '-w',
    # 'distance': '-d',
    'sv_fp': '-o',
    'clustal_fp': '',
}

value_config = {
    'tool_fp': 'python',
    'method': 'js_divergence',
    'window': '3',
    'distance': 'swissprot.distribution',
    'script_fpn': to('prot/feature/alignment/external/jsd/score_conservation.py'),
    'clustal_fp': to('data/msa/clustal/wild/SR24_AtoI/'),
    'sv_fp': to('data/jsd/SR24_AtoI/'),
}
```

In `prot.txt`, there are 7 proteins

``` text
ATAD2_LOC113841329
CAMK1G
CYP2W1_LOC101804267
KIF27
KIF27_LOC113841629
LOC119718710
RBBP8NL
```

Then, we can use `pp.external.jsd` to running JSD for a set of proteins.

``` py linenums="1"
import pypropel as pp

from pypropel.util.Reader import Reader as pfreader
df = pfreader().generic(df_fpn=to('data/msa/clustal/wild/SR24_AtoI/prot.txt'))
prots = df[0].unique()

for key, prot in enumerate(prots):
    order_list = [
        value_config['tool_fp'],
        value_config['script_fpn'],
        
        param_config['method'], value_config['method'],
        # param_config['window'], value_config['window'],
        # param_config['distance'], value_config['distance'],
        param_config['sv_fp'], value_config['sv_fp'] + prot + '.jsd',
        param_config['clustal_fp'], value_config['clustal_fp'] + prot + '.clustal',
    ]
    pp.external.jsd(
        order_list=order_list,
        job_fp='./',
        job_fn=str(key),
    )
```