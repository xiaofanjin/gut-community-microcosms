import pandas as pd
import numpy as np
import sys
import re
from glob import glob
import io
import os
from Bio import SearchIO
import gzip
import os
from functools import partial

def hmmOverlap(hit):
    sigDomains=[domain for domain in hit.hsps if domain.evalue<1e-2]
    hmmRanges=sorted([domain.query_range for domain in sigDomains])
    nonOverlappingRanges=[(0,0)]
    for ind in range(len(hmmRanges)):
        newRange=hmmRanges.pop(0)
        oldRange=nonOverlappingRanges[-1]
        if newRange[0]<=oldRange[1]:
            nonOverlappingRanges[-1]=(np.min([newRange[0],oldRange[0]]),np.max([newRange[1],oldRange[1]]))
        else:
            nonOverlappingRanges.append(newRange)
    return sum([interval[1]-interval[0] for interval in nonOverlappingRanges])
def kofamHmsToDf(hmsFile,outFile,makeHitsFiles=False):
    suffix=os.path.splitext(hmsFile)[-1]
    _open = partial(gzip.open, mode='rt') if suffix == '.gz' else open
    with _open(hmsFile) as f:
        contents = f.read()
    contentsSplit=re.split('^K[0-9]{5}',contents, flags=re.MULTILINE)
    data=[]
    for contentSplit in contentsSplit:
        for qresult in SearchIO.parse(io.StringIO(contentSplit), 'hmmer3-text'):
            pass
        try:
            if len(qresult)>0:
                data=data+[(hit.query_id,qresult.seq_len,hit.id,hit.description,hit.evalue,hit.bitscore,kolistDf.loc[kolistDf.knum==hit.query_id,'threshold'].to_numpy()[0],hmmOverlap(hit),hmmOverlap(hit)/qresult.seq_len) for hit in qresult.hits]
        except:
            pass
    df=pd.DataFrame(data, columns=['hmm_queryID','hmm_queryLen','hit_ID','hit_description','hit_evalue','hit_bitscore','bitscore_threshold','hit_overlap','hit_overlapFrac'])
    df['uhgg_ID']=os.path.basename(hmsFile).split('_')[0]
    df.to_csv(outFile,compression='gzip',index=False)
    if makeHitsFiles:
        bitscoreHits=df[df.bitscore_threshold<df.hit_bitscore].hmm_queryID.drop_duplicates()
        overlapHits=df[df.hit_overlapFrac>0.5].hmm_queryID.drop_duplicates()
        mixedHalfHits=df[(df.hit_overlapFrac>0.5) & (df.bitscore_threshold<2*df.hit_bitscore)].hmm_queryID.drop_duplicates()
        bitscoreHits.to_csv(outFile.split('.csv')[0]+'_bitscoreHits.txt.gz',index=False,header=False)
        overlapHits.to_csv(outFile.split('.csv')[0]+'_overlapHits.txt.gz',index=False,header=False)
        mixedHalfHits.to_csv(outFile.split('.csv')[0]+'_mixedHalfHits.txt.gz',index=False,header=False)
    return df

if __name__=="__main__":
    hmsFile=sys.argv[1]
    makeHitsFiles=False
    if len(sys.argv)>=3:
        outFile=sys.argv[2]
    else:
        outFile=hmsFile.split('.hms')[0]+'.csv.gz'
    if len(sys.argv)>=4:
        koListFile=sys.argv[3]
    else:
        koListFile='/pollard/data/protein_families/kegg/kofamProfiles20210204/ko_list' 
    if len(sys.argv)==5:
        makeHitsFiles=(sys.argv[4]=="makehits")
    kolistDf=pd.read_csv(koListFile,sep="\t")
    kolistDf['threshold']=pd.to_numeric(kolistDf['threshold'],errors='coerce')
    kofamHmsToDf(hmsFile,outFile,makeHitsFiles)
