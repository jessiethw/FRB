"""Setup FRB analysis
Jessie Thwaites 8/20/21
GFU, v2p6
"""

import numpy as np
import csky as cy

def load_new():
        import datetime
        date=datetime.datetime.now()

        ana_dir = cy.utils.ensure_dir('/home/jthwaites/csky_cache/'+date.strftime("%y-%m-%d"))
        ana = cy.get_analysis(cy.selections.repo, 'version-002-p06', cy.selections.GFUDataSpecs.GFU_IC86)

        conf = {'extended': True, #use extended LLH due to low time window
                'space': "ps",
                'time': "transient",
                'sig': 'transient',
                'ana':ana,
                'mp_cpus': 5, #some functions fail with >1 (?)
                'extra_keep':['energy']
                }

        cy.CONF.update(conf)

        ana.save(ana_dir)

def reload_ana():
        import glob
        import os

        dates=[os.path.basename(file) for file in glob.glob('/home/jthwaites/csky_cache/*')]

        if len(dates)==0:
                'No cached analysis object found, loading new'
                load_new()
                return
        
        most_recent=str(max([int(ymd.replace('-','')) for ymd in dates]))

        ana_dir = cy.utils.ensure_dir('/home/jthwaites/csky_cache/%s-%s-%s' \
                %(most_recent[:2],most_recent[2:4],most_recent[4:]))

        ana=cy.get_analysis(cy.selections.repo, 'version-002-p06', cy.selections.GFUDataSpecs.GFU_IC86, 
                        dir=ana_dir)

        conf = {'extended': True, #use extended LLH due to low time window
                'space': "ps",
                'time': "transient",
                'sig': 'transient',
                'ana':ana,
                'mp_cpus': 5, #some functions fail with >1 (?)
                'extra_keep':['energy']
                }

        cy.CONF.update(conf)
