"""Setup FRB analysis
Jessie Thwaites 8/20/21
GFU, v2p6
"""

import numpy as np
import csky as cy

ana_dir = cy.utils.ensure_dir('/home/jthwaites/csky_cache')
ana = cy.get_analysis(cy.selections.repo, 'version-002-p06', cy.selections.GFUDataSpecs.GFU_IC86, 
                dir=ana_dir)

conf = {'extended': True, #use extended LLH due to low time window
        'space': "ps",
        'time': "transient",
        'sig': 'transient',
        'ana':ana,
        'mp_cpus': 5 #some functions fail with >1 (?)
        }

cy.CONF.update(conf)