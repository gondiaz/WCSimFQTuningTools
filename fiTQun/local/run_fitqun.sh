#!/bin/bash

export FQTUNINGTOOLS=$HOME/HK_Software/WCSimFQTuningTools

. $HOME/HK_Software/fiTQun/install-Darwin_arm64-gcc_15.0.0-python_3.10.13/setup.sh

runfiTQun -p $FQTUNINGTOOLS/Time/WCTE_Parameters.dat -r out_fitqun.root $FQTUNINGTOOLS/Simulation/local/out/out.root
