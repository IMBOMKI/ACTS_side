#include <chrono>
#include <fstream>
#include <iostream>
#include <sstream>
#include <utility>
#include <assert.h>

#include <boost/type_erasure/any_cast.hpp>

#include "Acts/Seeding/BinFinder.hpp"
#include "Acts/Seeding/BinnedSPGroup.hpp"
#include "Acts/Seeding/InternalSeed.hpp"
#include "Acts/Seeding/InternalSpacePoint.hpp"
#include "Acts/Seeding/Seed.hpp"
#include "Acts/Seeding/SeedFilter.hpp"
#include "Acts/Seeding/Seedfinder.hpp"
#include "Acts/Seeding/SpacePointGrid.hpp"

#include "ATLASCuts.hpp"
#include "SpacePoint.hpp"
#include <Acts/Seeding/SeedfinderCPUFunctions.hpp>

// --- ROOT headers --- //
#include "TFile.h"
#include "TTree.h"

typedef SpacePoint external_spacepoint_t;
typedef Acts::Neighborhood<external_spacepoint_t> sp_range_t;
 
std::vector<const SpacePoint*> readFile(std::string filename) {
  std::string line;
  int layer;
  std::vector<const SpacePoint*> readSP;

  std::ifstream spFile(filename);
  if (spFile.is_open()) {
    while (!spFile.eof()) {
      std::getline(spFile, line);
      std::stringstream ss(line);
      std::string linetype;
      ss >> linetype;
      float x, y, z, r, varianceR, varianceZ;
      if (linetype == "lxyz") {
        ss >> layer >> x >> y >> z >> varianceR >> varianceZ;
        r = std::sqrt(x * x + y * y);
        float f22 = varianceR;
        float wid = varianceZ;
        float cov = wid * wid * .08333;
        if (cov < f22)
          cov = f22;
        if (std::abs(z) > 450.) {
          varianceZ = 9. * cov;
          varianceR = .06;
        } else {
          varianceR = 9. * cov;
          varianceZ = .06;
        }
        SpacePoint* sp =
            new SpacePoint{x, y, z, r, layer, varianceR, varianceZ};
        //     if(r < 200.){
        //       sp->setClusterList(1,0);
        //     }
        readSP.push_back(sp);
      }
    }
  }
  return readSP;
}

template <typename external_spacepoint_t>
void transformCoordinates(
			  std::vector<const Acts::InternalSpacePoint<external_spacepoint_t>*>& vec,
			  const Acts::InternalSpacePoint<external_spacepoint_t>& spM, bool bottom,
			  std::vector<Acts::LinCircle>& linCircleVec);

  
int main(int argc, char** argv) {
  //std::string file{"sp.txt"};
  std::string file{"sample_10k.txt"};
  bool help(false);
  bool quiet(false);

  int opt;
  while ((opt = getopt(argc, argv, "hf:q")) != -1) {
    switch (opt) {
      case 'f':
        file = optarg;
        break;
      case 'q':
        quiet = true;
        break;
      case 'h':
        help = true;
        [[fallthrough]];
      default: /* '?' */
        std::cerr << "Usage: " << argv[0] << " [-hq] [-f FILENAME]\n";
        if (help) {
          std::cout << "      -h : this help" << std::endl;
          std::cout
              << "      -f FILE : read spacepoints from FILE. Default is \""
              << file << "\"" << std::endl;
          std::cout << "      -q : don't print out all found seeds"
                    << std::endl;
        }

        exit(EXIT_FAILURE);
    }
  }

  std::ifstream f(file);
  if (!f.good()) {
    std::cerr << "input file \"" << file << "\" does not exist\n";
    exit(EXIT_FAILURE);
  }

  auto start_read = std::chrono::system_clock::now();
  std::vector<const SpacePoint*> spVec = readFile(file);
  auto end_read = std::chrono::system_clock::now();
  std::chrono::duration<double> elapsed_read = end_read - start_read;

  std::cout << "read " << spVec.size() << " SP from file " << file << " in "
            << elapsed_read.count() << "s" << std::endl;

  Acts::SeedfinderConfig<SpacePoint> config;
  // silicon detector max
  config.rMax = 160.;
  config.deltaRMin = 5.;
  config.deltaRMax = 160.;
  config.collisionRegionMin = -250.;
  config.collisionRegionMax = 250.;
  config.zMin = -2800.;
  config.zMax = 2800.;
  config.maxSeedsPerSpM = 5;
  // 2.7 eta
  config.cotThetaMax = 7.40627;
  config.sigmaScattering = 1.00000;

  config.minPt = 500.;
  config.bFieldInZ = 0.00199724;

  config.beamPos = {-.5, -.5};
  config.impactMax = 10.;

  auto bottomBinFinder = std::make_shared<Acts::BinFinder<SpacePoint>>(
      Acts::BinFinder<SpacePoint>());
  auto topBinFinder = std::make_shared<Acts::BinFinder<SpacePoint>>(
      Acts::BinFinder<SpacePoint>());
  Acts::SeedFilterConfig sfconf;
  Acts::ATLASCuts<SpacePoint> atlasCuts = Acts::ATLASCuts<SpacePoint>();
  config.seedFilter = std::make_unique<Acts::SeedFilter<SpacePoint>>(
      Acts::SeedFilter<SpacePoint>(sfconf, &atlasCuts));

  // covariance tool, sets covariances per spacepoint as required
  auto ct = [=](const SpacePoint& sp, float, float, float) -> Acts::Vector2D {
    return {sp.varianceR, sp.varianceZ};
  };
  
  // setup spacepoint grid config
  Acts::SpacePointGridConfig gridConf;
  gridConf.bFieldInZ = config.bFieldInZ;
  gridConf.minPt = config.minPt;
  gridConf.rMax = config.rMax;
  gridConf.zMax = config.zMax;
  gridConf.zMin = config.zMin;
  gridConf.deltaRMax = config.deltaRMax;
  gridConf.cotThetaMax = config.cotThetaMax;
  // create grid with bin sizes according to the configured geometry
  std::unique_ptr<Acts::SpacePointGrid<SpacePoint>> grid =
      Acts::SpacePointGridCreator::createGrid<SpacePoint>(gridConf);
  auto spGroup = Acts::BinnedSPGroup<SpacePoint>(spVec.begin(), spVec.end(), ct,
                                                 bottomBinFinder, topBinFinder,
                                                 std::move(grid), config);

  // w/ module
  
  std::vector<std::vector<Acts::Seed<SpacePoint>>> seedVector_wM;
  auto groupIt = spGroup.begin();
  auto endOfGroups = spGroup.end();

  for (; !(groupIt == endOfGroups); ++groupIt) {
    
    auto bottomSPs = groupIt.bottom();
    auto middleSPs = groupIt.middle();
    auto topSPs    = groupIt.top();

    std::vector<Acts::Seed<external_spacepoint_t>> outputVec;
    
    for (auto spM : middleSPs) {    
      
      float rM = spM->radius();
      float zM = spM->z();
      float varianceRM = spM->varianceR();
      float varianceZM = spM->varianceZ();
      
      // Doublet search    
      auto compatBottomSP =
	Acts::SeedfinderCPUFunctions<external_spacepoint_t,
				     sp_range_t>::SearchDoublet(true, bottomSPs, *spM, config);
      
      // no bottom SP found -> try next spM
      if (compatBottomSP.empty()) {
	continue;
      }
      
      auto compatTopSP =
	Acts::SeedfinderCPUFunctions<external_spacepoint_t,
				     sp_range_t>::SearchDoublet(false, topSPs, *spM, config);
      
      // no top SP found -> try next spM
      if (compatTopSP.empty()) {
	continue;
      }
      
      // contains parameters required to calculate circle with linear equation
      // ...for bottom-middle
      std::vector<Acts::LinCircle> linCircleBottom;
      // ...for middle-top
      std::vector<Acts::LinCircle> linCircleTop;
      
      Acts::SeedfinderCPUFunctions<external_spacepoint_t,sp_range_t>::transformCoordinates(compatBottomSP, *spM, true, linCircleBottom);
      Acts::SeedfinderCPUFunctions<external_spacepoint_t,sp_range_t>::transformCoordinates(compatTopSP, *spM, false, linCircleTop);
      
      auto seedsPerSpM = Acts::SeedfinderCPUFunctions<external_spacepoint_t,sp_range_t>::SearchTriplet(*spM, compatBottomSP, compatTopSP, linCircleBottom, linCircleTop, config);
      
      config.seedFilter->filterSeeds_1SpFixed(seedsPerSpM, outputVec);
      
    }

    seedVector_wM.push_back(outputVec);    
  }  


  // w/o module

  std::vector<std::vector<Acts::Seed<SpacePoint>>> seedVector_oM;

  groupIt = spGroup.begin();
  for (; !(groupIt == endOfGroups); ++groupIt) {
    
    auto bottomSPs = groupIt.bottom();
    auto middleSPs = groupIt.middle();
    auto topSPs    = groupIt.top();

    std::vector<Acts::Seed<external_spacepoint_t>> outputVec;

    /// Put Seedfinder.ipp contents in official repository
    for (auto spM : middleSPs) {
      float rM = spM->radius();
      float zM = spM->z();
      float varianceRM = spM->varianceR();
      float varianceZM = spM->varianceZ();
      
      // bottom space point
      std::vector<const Acts::InternalSpacePoint<external_spacepoint_t>*>
	compatBottomSP;
      
      for (auto bottomSP : bottomSPs) {
	float rB = bottomSP->radius();
	float deltaR = rM - rB;
	// if r-distance is too big, try next SP in bin
	if (deltaR > config.deltaRMax) {
	  continue;
	}
	// if r-distance is too small, break because bins are NOT r-sorted
	if (deltaR < config.deltaRMin) {
	  continue;
	}
	// ratio Z/R (forward angle) of space point duplet
	float cotTheta = (zM - bottomSP->z()) / deltaR;
	if (std::fabs(cotTheta) > config.cotThetaMax) {
	  continue;
	}
	// check if duplet origin on z axis within collision region
	float zOrigin = zM - rM * cotTheta;
	if (zOrigin < config.collisionRegionMin ||
	    zOrigin > config.collisionRegionMax) {
	  continue;
	}
	compatBottomSP.push_back(bottomSP);
      }
      // no bottom SP found -> try next spM
      if (compatBottomSP.empty()) {
	continue;
      }
      
      std::vector<const Acts::InternalSpacePoint<external_spacepoint_t>*> compatTopSP;
      
      for (auto topSP : topSPs) {
	float rT = topSP->radius();
	float deltaR = rT - rM;
	// this condition is the opposite of the condition for bottom SP
	if (deltaR < config.deltaRMin) {
	  continue;
	}
	if (deltaR > config.deltaRMax) {
	  break;
	}
	
	float cotTheta = (topSP->z() - zM) / deltaR;
	if (std::fabs(cotTheta) > config.cotThetaMax) {
	  continue;
	}
	float zOrigin = zM - rM * cotTheta;
	if (zOrigin < config.collisionRegionMin ||
	    zOrigin > config.collisionRegionMax) {
	  continue;
	  }
	compatTopSP.push_back(topSP);
      }
      if (compatTopSP.empty()) {
	continue;
      }
      // contains parameters required to calculate circle with linear equation
      // ...for bottom-middle
      std::vector<Acts::LinCircle> linCircleBottom;
      // ...for middle-top
      std::vector<Acts::LinCircle> linCircleTop;
      transformCoordinates<external_spacepoint_t>(compatBottomSP, *spM, true, linCircleBottom);
      transformCoordinates<external_spacepoint_t>(compatTopSP, *spM, false, linCircleTop);
      
      // create vectors here to avoid reallocation in each loop
      std::vector<const Acts::InternalSpacePoint<external_spacepoint_t>*> topSpVec;
      std::vector<float> curvatures;
      std::vector<float> impactParameters;
      
      std::vector<std::pair<
	float, std::unique_ptr<const Acts::InternalSeed<external_spacepoint_t>>>>
	seedsPerSpM;
      size_t numBotSP = compatBottomSP.size();
      size_t numTopSP = compatTopSP.size();
	
      for (size_t b = 0; b < numBotSP; b++) {
	auto lb = linCircleBottom[b];
	float Zob = lb.Zo;
	float cotThetaB = lb.cotTheta;
	float Vb = lb.V;
	float Ub = lb.U;
	float ErB = lb.Er;
	float iDeltaRB = lb.iDeltaR;
	
	// 1+(cot^2(theta)) = 1/sin^2(theta)
	float iSinTheta2 = (1. + cotThetaB * cotThetaB);
	// calculate max scattering for min momentum at the seed's theta angle
	// scaling scatteringAngle^2 by sin^2(theta) to convert pT^2 to p^2
	// accurate would be taking 1/atan(thetaBottom)-1/atan(thetaTop) <
	// scattering
	// but to avoid trig functions we approximate cot by scaling by
	// 1/sin^4(theta)
	// resolving with pT to p scaling --> only divide by sin^2(theta)
	// max approximation error for allowed scattering angles of 0.04 rad at
	// eta=infinity: ~8.5%
	float scatteringInRegion2 = config.maxScatteringAngle2 * iSinTheta2;
	// multiply the squared sigma onto the squared scattering
	scatteringInRegion2 *=
	  config.sigmaScattering * config.sigmaScattering;
	
	// clear all vectors used in each inner for loop
	topSpVec.clear();
	curvatures.clear();
	impactParameters.clear();
	for (size_t t = 0; t < numTopSP; t++) {
	  auto lt = linCircleTop[t];
	  
	  // add errors of spB-spM and spM-spT pairs and add the correlation term
	  // for errors on spM
	  float error2 = lt.Er + ErB +
	    2 * (cotThetaB * lt.cotTheta * varianceRM + varianceZM) *
	    iDeltaRB * lt.iDeltaR;
	  
	  float deltaCotTheta = cotThetaB - lt.cotTheta;
	  float deltaCotTheta2 = deltaCotTheta * deltaCotTheta;
	  float error;
	  float dCotThetaMinusError2;
	  // if the error is larger than the difference in theta, no need to
	  // compare with scattering
	  if (deltaCotTheta2 - error2 > 0) {
	    deltaCotTheta = std::abs(deltaCotTheta);
	    // if deltaTheta larger than the scattering for the lower pT cut, skip
	    error = std::sqrt(error2);
	    dCotThetaMinusError2 =
	      deltaCotTheta2 + error2 - 2 * deltaCotTheta * error;
	    // avoid taking root of scatteringInRegion
	    // if left side of ">" is positive, both sides of unequality can be
	    // squared
	    // (scattering is always positive)
	    
	    if (dCotThetaMinusError2 > scatteringInRegion2) {
	      continue;
	    }
	  }
	    
	  // protects against division by 0
	  float dU = lt.U - Ub;
	  if (dU == 0.) {
	    continue;
	  }
	  // A and B are evaluated as a function of the circumference parameters
	  // x_0 and y_0
	  float A = (lt.V - Vb) / dU;
	  float S2 = 1. + A * A;
	  float B = Vb - A * Ub;
	  float B2 = B * B;
	  // sqrt(S2)/B = 2 * helixradius
	  // calculated radius must not be smaller than minimum radius
	  if (S2 < B2 * config.minHelixDiameter2) {
	    continue;
	  }
	  // 1/helixradius: (B/sqrt(S2))/2 (we leave everything squared)
	  float iHelixDiameter2 = B2 / S2;
	  // calculate scattering for p(T) calculated from seed curvature
	  float pT2scatter = 4 * iHelixDiameter2 * config.pT2perRadius;
	  // TODO: include upper pT limit for scatter calc
	  // convert p(T) to p scaling by sin^2(theta) AND scale by 1/sin^4(theta)
	  // from rad to deltaCotTheta
	  float p2scatter = pT2scatter * iSinTheta2;
	  // if deltaTheta larger than allowed scattering for calculated pT, skip
	  if ((deltaCotTheta2 - error2 > 0) &&
	      (dCotThetaMinusError2 >
	       p2scatter * config.sigmaScattering * config.sigmaScattering)) {
	    continue;
	  }
	  // A and B allow calculation of impact params in U/V plane with linear
	  // function
	  // (in contrast to having to solve a quadratic function in x/y plane)
	  float Im = std::abs((A - B * rM) * rM);
	  
	  if (Im <= config.impactMax) {
	    topSpVec.push_back(compatTopSP[t]);
	    // inverse diameter is signed depending if the curvature is
	    // positive/negative in phi
	    curvatures.push_back(B / std::sqrt(S2));
	    impactParameters.push_back(Im);
	  }
	}
	if (!topSpVec.empty()) {
	  std::vector<std::pair<
	    float, std::unique_ptr<const Acts::InternalSeed<external_spacepoint_t>>>>
	    sameTrackSeeds;
	  sameTrackSeeds = std::move(config.seedFilter->filterSeeds_2SpFixed(
									       *compatBottomSP[b], *spM, topSpVec, curvatures, impactParameters,
									       Zob));
	  seedsPerSpM.insert(seedsPerSpM.end(),
			     std::make_move_iterator(sameTrackSeeds.begin()),
			     std::make_move_iterator(sameTrackSeeds.end()));
	}
      }
      config.seedFilter->filterSeeds_1SpFixed(seedsPerSpM, outputVec);
    }
    
     
    seedVector_oM.push_back(outputVec);    
  }  

  // Assert the sizes from two methods are same
  assert(seedVector_wM.size() == seedVector_oM.size());
  
  int numSeeds = 0;
  for (auto& outVec : seedVector_wM) {
    numSeeds += outVec.size();
  }
  std::cout << "Number of seeds generated: " << numSeeds << std::endl;
  if (!quiet) {
    
    for (size_t i=0; i<numSeeds; i++){
      
      auto& regionVec_wM = seedVector_wM[i];
      auto& regionVec_oM = seedVector_oM[i];

      std::cout << regionVec_wM.size() << "  " << regionVec_oM.size() << std::endl;
      assert(regionVec_wM.size() == regionVec_oM.size());
      
      for (size_t j=0; j<regionVec_wM.size(); j++){
	
	const Acts::Seed<SpacePoint>* seed_wM = &regionVec_wM[j];
        const SpacePoint* sp_wM = seed_wM->sp()[0];

	const Acts::Seed<SpacePoint>* seed_oM = &regionVec_oM[j];
        const SpacePoint* sp_oM = seed_oM->sp()[0];

	assert(sp_wM->x() == sp_oM->x());
	assert(sp_wM->y() == sp_oM->y());
	assert(sp_wM->z() == sp_oM->z());
	
	std::cout << " (" << sp_wM->x() << ", " << sp_wM->y() << ", " << sp_wM->z()
                  << ") ";
        sp_wM = seed_wM->sp()[1];
        std::cout << sp_wM->surface << " (" << sp_wM->x() << ", " << sp_wM->y() << ", "
                  << sp_wM->z() << ") ";
        sp_wM = seed_wM->sp()[2];
        std::cout << sp_wM->surface << " (" << sp_wM->x() << ", " << sp_wM->y() << ", "
                  << sp_wM->z() << ") ";
        std::cout << std::endl;
	
      } 
      
    }
    
    /*
    for (auto& regionVec : seedVector_wM) {
      for (size_t i = 0; i < regionVec.size(); i++) {
        const Acts::Seed<SpacePoint>* seed = &regionVec[i];
        const SpacePoint* sp = seed->sp()[0];
        std::cout << " (" << sp->x() << ", " << sp->y() << ", " << sp->z()
                  << ") ";
        sp = seed->sp()[1];
        std::cout << sp->surface << " (" << sp->x() << ", " << sp->y() << ", "
                  << sp->z() << ") ";
        sp = seed->sp()[2];
        std::cout << sp->surface << " (" << sp->x() << ", " << sp->y() << ", "
                  << sp->z() << ") ";
        std::cout << std::endl;
      }
    }
    */
  }
  
  return 0;
}
  


template <typename external_spacepoint_t>
void transformCoordinates(
    std::vector<const Acts::InternalSpacePoint<external_spacepoint_t>*>& vec,
    const Acts::InternalSpacePoint<external_spacepoint_t>& spM, bool bottom,
    std::vector<Acts::LinCircle>& linCircleVec) {
  float xM = spM.x();
  float yM = spM.y();
  float zM = spM.z();
  float rM = spM.radius();
  float varianceZM = spM.varianceZ();
  float varianceRM = spM.varianceR();
  float cosPhiM = xM / rM;
  float sinPhiM = yM / rM;
  for (auto sp : vec) {
    float deltaX = sp->x() - xM;
    float deltaY = sp->y() - yM;
    float deltaZ = sp->z() - zM;
    // calculate projection fraction of spM->sp vector pointing in same
    // direction as
    // vector origin->spM (x) and projection fraction of spM->sp vector pointing
    // orthogonal to origin->spM (y)
    float x = deltaX * cosPhiM + deltaY * sinPhiM;
    float y = deltaY * cosPhiM - deltaX * sinPhiM;
    // 1/(length of M -> SP)
    float iDeltaR2 = 1. / (deltaX * deltaX + deltaY * deltaY);
    float iDeltaR = std::sqrt(iDeltaR2);
    //
    int bottomFactor = 1 * (int(!bottom)) - 1 * (int(bottom));
    // cot_theta = (deltaZ/deltaR)
    float cot_theta = deltaZ * iDeltaR * bottomFactor;
    // VERY frequent (SP^3) access
    Acts::LinCircle l;
    l.cotTheta = cot_theta;
    // location on z-axis of this SP-duplet
    l.Zo = zM - rM * cot_theta;
    l.iDeltaR = iDeltaR;
    // transformation of circle equation (x,y) into linear equation (u,v)
    // x^2 + y^2 - 2x_0*x - 2y_0*y = 0
    // is transformed into
    // 1 - 2x_0*u - 2y_0*v = 0
    // using the following m_U and m_V
    // (u = A + B*v); A and B are created later on
    l.U = x * iDeltaR2;
    l.V = y * iDeltaR2;
    // error term for sp-pair without correlation of middle space point
    l.Er = ((varianceZM + sp->varianceZ()) +
            (cot_theta * cot_theta) * (varianceRM + sp->varianceR())) *
           iDeltaR2;
    linCircleVec.push_back(l);
  }
}
