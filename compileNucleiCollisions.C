void compileNucleiCollisions()
{
    gROOT->ProcessLine(".L StringGeneration/DistanceEntry.cpp+");
    gROOT->ProcessLine(".L StringGeneration/MinDistanceFinder.cpp+");
    gROOT->ProcessLine(".L StringGeneration/NucleiCollision.cpp+");

    gROOT->ProcessLine(".L ManagerNucleiCollisions.cpp+");

    gROOT->ProcessLine(".q");

}
