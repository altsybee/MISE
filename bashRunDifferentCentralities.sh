#!/bin/bash


echo "#####  start..."

rm points_*
rm -rf tmpOutputs
#rm canv*
#rm statistics*
#rm outputAnalysers_*
#int what_to_draw, int winId, int centralityStep, int listStart, int listEnd, const char * strOutName )

#inputFile='/Users/macbook/alice/aliceAnalysis/results/task_2014_07_13_PbPb_etaWins_Data_AOD_centrScan_blocks23_noHist2D_eW04_phi8/MergedOutput.root'

#tmp useful stuff:
#a=3
#b[0]=0
#b[1]=10
#b[2]=20
#for i in {0..2}
#        echo "#####  tests... $i... ${b[$i]}"



nEvents=10000

#for b in {0,5,10,15,20,25,30,35,40,45,50,55,60,65,70,75,80,85,90}
#for clusterFormDist in {0.1,0.15,0.2,0.25,0.3,0.35,0.4,0.45,0.5}

#for nPartons in {500,750,1000,1250,1500,1750,2000,2250,2500,2750,3000}
#for clusterFormDist in {0.55,0.6,0.65,0.7}
#for clusterFormDist in {0.1,0.15,0.2,0.25,0.3,0.35,0.4,0.45,0.5}

#for partonInteractionDist in {0.1,0.15,0.2,0.25,0.3,0.35,0.4,0.45,0.5,0.55,0.6,0.65,0.7}


#for meanPartonsInNucleon in {10,12,15} #5,8,18}
for meanPartonsInNucleon in {8,} #15 10
#for stringInteractionDist in {0.25,1,6}
#for coeffPtKickPerUnitMagn in {0.02,0.1,0.2}
do
    partonInteractionDist=0.25
    clusterFormDist=20 #0.4 #0.3
    stringInteractionDist=2 #6 #3
    strIntDistFractionParA=5.
    coeffPtKickPerUnitMagn=0.25 #0.05
#    meanPartonsInNucleon=15

    fileLabel=pDist${partonInteractionDist}_strDist${stringInteractionDist}_clust${clusterFormDist}_nP${meanPartonsInNucleon}_kick${coeffPtKickPerUnitMagn}
    mkdir out_${fileLabel}

    #for b in {4,5} #for tests
    #for b in {0,6,11} #1,2,3,4,5,6,7,8,9,10,11,12,13,14,15}
    for b in {0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15}
    do
        echo "#####  impact parameter $b"
        root -l "launchNuclearStructure.C( ${nEvents}, $b, $partonInteractionDist, $clusterFormDist, $stringInteractionDist, $strIntDistFractionParA, $meanPartonsInNucleon, $coeffPtKickPerUnitMagn, 2 )"

        cat "tmpOutputs/tmpTextOutput0.txt" >> points_percolationParameter_${fileLabel}.txt
        cat "tmpOutputs/tmpTextOutput1.txt" >> points_nStrRatio_${fileLabel}.txt
        cat "tmpOutputs/tmpTextOutput1_0.txt" >> points_nStrings_${fileLabel}.txt
        cat "tmpOutputs/tmpTextOutput2.txt" >> points_Npart_${fileLabel}.txt
        cat "tmpOutputs/tmpTextOutput3.txt" >> points_Ncoll_${fileLabel}.txt
        cat "tmpOutputs/tmpTextOutput4.txt" >> points_FxFyAssym_${fileLabel}.txt
        cat "tmpOutputs/tmpTextOutput5.txt" >> points_strInteractionMagn_${fileLabel}.txt

        cat "tmpOutputs/tmpTextOutput_EventManager_ptVsNch.txt" >> points_meanPt_${fileLabel}.txt


  cat "tmpOutputs/tmpTextOutput_bNN_list0_eta0_phi0.txt" >> points_bNN_list0_eta0_phi0_${fileLabel}.txt
  cat "tmpOutputs/tmpTextOutput_bNN_list0_eta0_phi1.txt" >> points_bNN_list0_eta0_phi1_${fileLabel}.txt
  cat "tmpOutputs/tmpTextOutput_bNN_list0_eta0_phi2.txt" >> points_bNN_list0_eta0_phi2_${fileLabel}.txt
  cat "tmpOutputs/tmpTextOutput_bNN_list0_eta0_phi3.txt" >> points_bNN_list0_eta0_phi3_${fileLabel}.txt
  cat "tmpOutputs/tmpTextOutput_bNN_list0_eta0_phi4.txt" >> points_bNN_list0_eta0_phi4_${fileLabel}.txt
  cat "tmpOutputs/tmpTextOutput_bNN_list0_eta0_phi5.txt" >> points_bNN_list0_eta0_phi5_${fileLabel}.txt
  cat "tmpOutputs/tmpTextOutput_bNN_list0_eta0_phi6.txt" >> points_bNN_list0_eta0_phi6_${fileLabel}.txt
  cat "tmpOutputs/tmpTextOutput_bNN_list0_eta0_phi7.txt" >> points_bNN_list0_eta0_phi7_${fileLabel}.txt
  cat "tmpOutputs/tmpTextOutput_bNN_list1_eta0_phi0.txt" >> points_bNN_list1_eta0_phi0_${fileLabel}.txt
  cat "tmpOutputs/tmpTextOutput_bNN_list1_eta0_phi1.txt" >> points_bNN_list1_eta0_phi1_${fileLabel}.txt
  cat "tmpOutputs/tmpTextOutput_bNN_list1_eta0_phi2.txt" >> points_bNN_list1_eta0_phi2_${fileLabel}.txt
  cat "tmpOutputs/tmpTextOutput_bNN_list1_eta0_phi3.txt" >> points_bNN_list1_eta0_phi3_${fileLabel}.txt
  cat "tmpOutputs/tmpTextOutput_bNN_list1_eta0_phi4.txt" >> points_bNN_list1_eta0_phi4_${fileLabel}.txt
  cat "tmpOutputs/tmpTextOutput_bNN_list1_eta0_phi5.txt" >> points_bNN_list1_eta0_phi5_${fileLabel}.txt
  cat "tmpOutputs/tmpTextOutput_bNN_list1_eta0_phi6.txt" >> points_bNN_list1_eta0_phi6_${fileLabel}.txt
  cat "tmpOutputs/tmpTextOutput_bNN_list1_eta0_phi7.txt" >> points_bNN_list1_eta0_phi7_${fileLabel}.txt
  cat "tmpOutputs/tmpTextOutput_bNN_list2_eta0_phi0.txt" >> points_bNN_list2_eta0_phi0_${fileLabel}.txt
  cat "tmpOutputs/tmpTextOutput_bNN_list3_eta0_phi0.txt" >> points_bNN_list3_eta0_phi0_${fileLabel}.txt
cat "tmpOutputs/tmpTextOutput_bPtPt_list0_eta0_phi0.txt" >> points_bPtPt_list0_eta0_phi0_${fileLabel}.txt
cat "tmpOutputs/tmpTextOutput_bPtPt_list0_eta0_phi1.txt" >> points_bPtPt_list0_eta0_phi1_${fileLabel}.txt
cat "tmpOutputs/tmpTextOutput_bPtPt_list0_eta0_phi2.txt" >> points_bPtPt_list0_eta0_phi2_${fileLabel}.txt
cat "tmpOutputs/tmpTextOutput_bPtPt_list0_eta0_phi3.txt" >> points_bPtPt_list0_eta0_phi3_${fileLabel}.txt
cat "tmpOutputs/tmpTextOutput_bPtPt_list0_eta0_phi4.txt" >> points_bPtPt_list0_eta0_phi4_${fileLabel}.txt
cat "tmpOutputs/tmpTextOutput_bPtPt_list0_eta0_phi5.txt" >> points_bPtPt_list0_eta0_phi5_${fileLabel}.txt
cat "tmpOutputs/tmpTextOutput_bPtPt_list0_eta0_phi6.txt" >> points_bPtPt_list0_eta0_phi6_${fileLabel}.txt
cat "tmpOutputs/tmpTextOutput_bPtPt_list0_eta0_phi7.txt" >> points_bPtPt_list0_eta0_phi7_${fileLabel}.txt
cat "tmpOutputs/tmpTextOutput_bPtPt_list1_eta0_phi0.txt" >> points_bPtPt_list1_eta0_phi0_${fileLabel}.txt
cat "tmpOutputs/tmpTextOutput_bPtPt_list1_eta0_phi1.txt" >> points_bPtPt_list1_eta0_phi1_${fileLabel}.txt
cat "tmpOutputs/tmpTextOutput_bPtPt_list1_eta0_phi2.txt" >> points_bPtPt_list1_eta0_phi2_${fileLabel}.txt
cat "tmpOutputs/tmpTextOutput_bPtPt_list1_eta0_phi3.txt" >> points_bPtPt_list1_eta0_phi3_${fileLabel}.txt
cat "tmpOutputs/tmpTextOutput_bPtPt_list1_eta0_phi4.txt" >> points_bPtPt_list1_eta0_phi4_${fileLabel}.txt
cat "tmpOutputs/tmpTextOutput_bPtPt_list1_eta0_phi5.txt" >> points_bPtPt_list1_eta0_phi5_${fileLabel}.txt
cat "tmpOutputs/tmpTextOutput_bPtPt_list1_eta0_phi6.txt" >> points_bPtPt_list1_eta0_phi6_${fileLabel}.txt
cat "tmpOutputs/tmpTextOutput_bPtPt_list1_eta0_phi7.txt" >> points_bPtPt_list1_eta0_phi7_${fileLabel}.txt
cat "tmpOutputs/tmpTextOutput_bPtPt_list2_eta0_phi0.txt" >> points_bPtPt_list2_eta0_phi0_${fileLabel}.txt
cat "tmpOutputs/tmpTextOutput_bPtPt_list3_eta0_phi0.txt" >> points_bPtPt_list3_eta0_phi0_${fileLabel}.txt


        # move histos and canvases
        mkdir out_${fileLabel}/b_${b}/
        mv tmpOutputs/* out_${fileLabel}/b_${b}/
        #mv canv*.* out_${fileLabel}/b_${b}/
        #mv statistics*.* out_${fileLabel}/b_${b}/
        #mv outputAnalysers* out_${fileLabel}/b_${b}/

    done
    # move data files
    mv points* out_${fileLabel}/

    echo "$nEvents" > nEvents.txt
    mv nEvents.txt out_${fileLabel}/

done
#root -l "launchNuclearStructure.C( $winId )"

rm -rf tmpOutputs


    #for b in {0,10,20,30,40,50,60,70,80,90}
    #for b in {0,5,10,15,20,25,30,35,40,45,50,55,60,65,70,75,80,85,90}

