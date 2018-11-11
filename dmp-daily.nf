#!/usr/bin/env nextflow

logPath = '$HOME/log'

/* Daily pull of MSK-IMPACT repository */
process pullRepo {
  
  script:
  """
    echo '\n' >> $logPath/msk-impact.log && \
    date -u >> $logPath/msk-impact.log && \
    hg pull -u -R /home/jonssonp/res/dmp &>> $logPath/msk-impact.log
  """
}

// process pairBamKey {

// }

// process queryGlioma {

// }

// process queryBRCA {

// }

// process newFacets {

// }

// process newSignatures {

// }