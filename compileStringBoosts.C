void compileStringBoosts()
{
    gROOT->ProcessLine(".L StringGeneration/DistanceEntry.cpp+");
    gROOT->ProcessLine(".L StringGeneration/MinDistanceFinder.cpp+");
    gROOT->ProcessLine(".L StringGeneration/StringBoosting.cpp+");


    gROOT->ProcessLine(".L ManagerStringBoostsCalculator.cpp+");

    gROOT->ProcessLine(".q");

}
