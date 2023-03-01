#include <cmath>
#include <fstream>
#include <iomanip>
#include <iostream>
#include "cuts.h"


bool ICenergyCut(
        SimpleDST dst, 
        photospline::splinetable<> &spline, 
        double zenith, double emin, double emax) {

  // Setup basic parameters
  double x = cos(zenith);
  double y = log10(dst.NChannels);

  // Boundary check (energy cut tables go to 0.3 in cos(zenith))
  if (x < 0.3)
    return false;

  // Catch additional outliers
  double coords[2] = {x, y};
  int centers[spline.get_ndim()];
  if (!spline.searchcenters(coords, centers)) {
    std::cout << "Variables outside of table boundaries" << std::endl;
    std::cout << "x: " << x << " y: " << y << std::endl;
    return false;
  }

  // Calculate reconstructed energy
  double median = spline.ndsplineeval(coords, centers, 0);
  // Make sure we're in the energy bin range
  if (emax <= 0.0)
     return (median >= emin); 
  return ((median >= emin) && (median < emax));
}

// Vector form
int ICenergyCut(SimpleDST dst, photospline::splinetable<> &spline, double zenith, std::vector<float> ebins) {

  // Setup basic parameters
  double x = cos(zenith);
  double y = log10(dst.NChannels);

  // Boundary check (energy cut tables go to 0.3 in cos(zenith))
  if (x < 0.3)
    return -1;

  // Catch additional outliers
  double coords[2] = {x, y};
  int centers[spline.get_ndim()];
  if (!spline.searchcenters(coords, centers)) {
    std::cout << "Variables outside of table boundaries" << std::endl;
    std::cout << "x: " << x << " y: " << y << std::endl;
    return -1;
  }

  // Calculate reconstructed energy
  double median = spline.ndsplineeval(coords, centers, 0);
  // Make sure we're in the energy bin range
  if ((median < ebins[0]) || (median > ebins.back()))
    return -1;

  // Get energy bin
  int ebin = 0;
  while (median > ebins[ebin+1])
    ebin += 1;

  return ebin;
}


bool ITenergyCut(SimpleDST dst, double emin, double emax) {

  // Get most likely energy value
  double llhEnergy = (dst.pLLH >= dst.fLLH) ? dst.pEnergy : dst.fEnergy;
  double logEnergy = log10(llhEnergy);

  if (emax <= 0.0)
    return (logEnergy >= emin);
  // Get energy bin
  return ((logEnergy >= emin) && (logEnergy < emax));
}

//vector form
int ITenergyCut(SimpleDST dst, std::vector<float> ebins) {

  // Get most likely energy value
  double llhEnergy = (dst.pLLH >= dst.fLLH) ? dst.pEnergy : dst.fEnergy;
  double logEnergy = log10(llhEnergy);

  // Make sure we're in the energy bin range
  if ((logEnergy < ebins[0]) || (logEnergy > ebins.back()))
    return -1;

  // Get energy bin
  int ebin = 0;
  while (logEnergy > ebins[ebin+1])
    ebin += 1;

  return ebin;
}


bool ITs125Cut(SimpleDST dst, double smin, double smax) {

  // Get desired s125 value
  double s125 = (dst.nStations >= 5) ? dst.s125 : dst.ss125;
  double logS125 = log10(s125);

  // Make sure we're in the bin range
  if (smax <= 0.0)
    return (logS125 >= smin);
  return ((logS125 >= smin) && (logS125 < smax));

}

// vector form
int ITs125Cut(SimpleDST dst, std::vector<float> sbins) {

  // Get desired s125 value
  double s125 = (dst.nStations >= 5) ? dst.s125 : dst.ss125;
  double logS125 = log10(s125);

  // Make sure we're in the bin range
  if ((logS125 < sbins[0]) || (logS125 > sbins.back()))
    return -1;

  // Get s125 bin
  int sbin = 0;
  while (logS125 > sbins[sbin+1])
    sbin += 1;

  return sbin;
}


bool ITNstatCut(SimpleDST dst, double smin, double smax) {
  //high energy is 8<numStat
  if (!(smax > 0))
    return (dst.nStations >=smin);
    
  return ( (dst.nStations >=smin) && (dst.nStations < smax));

}

int ITNstatCut(SimpleDST dst, std::vector<float> sbins) {
  if ((dst.nStations < sbins[0]) || (dst.nStations> sbins.back()))
    return -1;

  int sbin = 0;
  while (dst.nStations > sbins[sbin+1])
    sbin += 1;

  return sbin;

}


Config::Config(Detector dt, Configuration c, Filter f, Method m):
  detector(dt),
  cfg(c),
  filter(f) {}


Config::Config(po::variables_map vm)
{
  std::string config;
  if (vm.count("config") )
    config = vm["config"].as<std::string>();
  
  if (config=="IT59" || config=="IT73" || config=="IT81") {
    std::cout << "IceTop v1" << std::endl;
    detector = IceTop;
    cfg = ITv1;
  }
  else if (config=="IT81-2012" || config=="IT81-2013") {
    std::cout << "IceTop v2" << std::endl;
    detector = IceTop;
    cfg = ITv2;
  } 
  else if (config=="ITpass2") {
    std::cout << "IceTop v3" << std::endl;
    detector = IceTop;
    cfg = ITv3;
  } 
  else if (config=="IC86-2011" || config=="IC86-2012" || config=="IC86-2013" || 
      config=="IC86-2014" || config=="IC86-2015") 
  {
    std::cout << "IceCube v1" << std::endl;
    detector = IceCube;
    cfg = ICv1;
  } 
  else {
    std::cout << "IceCube v2" << std::endl;
    detector = IceCube;
    cfg = ICv2;
  } 

  std::string f;
  if (vm.count("filter") )
      f = vm["filter"].as<std::string>();
  if (f == "STA3")
      filter = STA3;
  else if (f == "STA8")
      filter = STA8;
  else if (f == "NotSTA8")
      filter = NotSTA8;
  else 
      filter = None;

  std::string m;
  if (vm.count("method") )
    m = vm["method"].as<std::string>();
  if (m == "ext") {
      method = extsid;
  }
  else if (m == "anti") {
      method = antisid;
  }
  else if (m == "solar") {
      method = solar;
  }
  else {
      method = sidereal;
  }
}


bool Config::newConfig() {
  return (cfg == ICv2);
}

bool filterCut(Config cfg, SimpleDST dst) {

  if (cfg.filter==Config::None)
     return true; 
  if (cfg.filter==Config::STA3 && dst.isSTA3)
    return true;
  if (cfg.cfg == Config::ITv1) {
    if ((cfg.filter==Config::STA8 && dst.isSTA8) ||
        (cfg.filter==Config::NotSTA8 && dst.isSTA3 && !dst.isSTA8))
      return true;
  }
  if (cfg.cfg==Config::ITv2) {
    if ((cfg.filter==Config::STA8 && dst.nStations>=8) ||
        (cfg.filter==Config::NotSTA8 && dst.isSTA3 && dst.nStations<8))
      return true;
  }
  if (cfg.cfg==Config::ITv3) {
    if ((cfg.filter==Config::STA8 && dst.nStations>=8) ||
        (cfg.filter==Config::NotSTA8 && dst.isSTA3 && dst.nStations<8))
      return true;
  }
  return false;
}

