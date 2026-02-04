void compileStringFragmentation()
{
    gROOT->ProcessLine(".L StringDecayer/ParticleDescr.cpp+");
    gROOT->ProcessLine(".L StringDecayer/DecayInTwo.cpp+");
    gROOT->ProcessLine(".L StringDecayer/StringFragmentation.cpp+");
    gROOT->ProcessLine(".L StringDecayer/StringDescr.cpp+");

    gROOT->ProcessLine(".L ManagerStringFragmentation.cpp+");

//    gROOT->ProcessLine(".q");

}
