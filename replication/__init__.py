"""
Python replication of Hasbrouck, "Price Discovery in High Resolution".

Shared code for all analyses lives in this package:
    helpers.py       trading-day iteration, chunked pulls, run listing (log) files
    taq_queries.py   WRDS TAQ (millisecond) queries
    nbbo.py          NBBO construction from exchange quotes (taqExtract01.sas) and getPrice.m
    mvar.py          sparse PDL VAR/VECM estimation and impulse responses (MVARi.m)
    random_walk.py   random-walk decomposition and information shares (randomWalkDecomp.m)
    wrds_batch.py    submit/monitor/fetch WRDS Cloud batch jobs from your machine
    wrds_job.sh      generic Grid Engine wrapper used by every analysis

Analysis-specific code lives in replication/analysisN/.
"""
