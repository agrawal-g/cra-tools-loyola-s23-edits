#include <LaputopSmall.h>
#include <string>
#include <TBranch.h>
#include <TChain.h>

void
LaputopSmall::SetupChain(TChain* chain, std::string config)
{
  std::string detector = config.substr(0,2);

  if (detector == "ITpass2") {

    chain->SetBranchAddress("time", &time, &b_time);
    chain->SetBranchAddress("Zenith", &zenith, &b_zenith);
    chain->SetBranchAddress("Azimuth", &azimuth, &b_azimuth);
  }

}
