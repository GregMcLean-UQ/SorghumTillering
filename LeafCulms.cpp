//------------------------------------------------------------------------------------------------
#include <stdio.h>
#include "Plant.h"
#include "LeafCulms.h"

using namespace Sorghum;

//------------------------------------------------------------------------------------------------
//------ LeafCulms
//------------------------------------------------------------------------------------------------

//------------------------------------------------------------------------------------------------
//------ LeafCulms Constructor
//------------------------------------------------------------------------------------------------
LeafCulms::LeafCulms(ScienceAPI2& api, Plant* p) : Leaf(api, p)
{
	plant = p;
	name = "Leaf";
	partNo = 1;
	//verticalAdjustment = 0.1;
	//aMaxVert = 0.3;
	//aTillerVert = 0.05;
	avgRadiation = 0.0;
	calculatedTillers = 0.0;
	thermalTimeCount = 0.0;
	maxLAIForTillerAddition = 0.325;
	startThermalQuotientLeafNo = 3;
	endThermalQuotientLeafNo = 5;
	linearLAI = 0.0;
	initialize();
	doRegistrations();

}
//------------------------------------------------------------------------------------------------
//------ Destructor
//------------------------------------------------------------------------------------------------
LeafCulms::~LeafCulms()
{
	//delete culms
	while (!Culms.empty())
		delete Culms.back(), Culms.pop_back();
}

void LeafCulms::initialize()
{
	//Initialize shouldn't be called once sown.
	//Remove all then add the first culm (which is the main culm or stem).

	Culms.clear();

	// Initialise Main Stem
	initiateTiller(0, 1, 1);

	tillersAdded = 0;
	calculatedTillers = 0.0;
	supply = 0;
	demand = 0;
	tillers = 0.0;
	tillerSdIntercept = 0;
	tillerSdSlope = 0;
	tillerSlaBound = 0;
	radiationValues = 0;
	temperatureValues = 0;
	laiReductionForSLA = 0.0;
	totalLaiReductionForSLA = 0.0;
	maxLaiTarget = 0.0;
	tillerLaiToReduce = 0.0;

	Leaf::initialize();
}
//------------------------------------------------------------------------------------------------
//-----------  read leaf parameters
//------------------------------------------------------------------------------------------------
void LeafCulms::readParams(void)
{
	Leaf::readParams();
	// leaf area individual leaf
	//Birch, Hammer bell shaped curve parameters
	scienceAPI.read("aX0S", "", false, leafAreaParams.aX0S);    // Eqn 14
	scienceAPI.read("aX0I", "", false, leafAreaParams.aX0I);    // Eqn 14 adj EVO Apr 2023
	scienceAPI.read("aMaxI", "", false, leafAreaParams.aMaxI); // Eqn 13 adj EVO Apr 2023
	scienceAPI.read("aMaxS", "", false, leafAreaParams.aMaxS);
	scienceAPI.read("leaf_no_correction", "", false, leafAreaParams.leafNoCorrection);
	scienceAPI.read("largestLeafPlateau", "", true, leafAreaParams.largestLeafPlateau);

	// DEBUG Make External
	leafAreaParams.a0 = 0.00955;
	leafAreaParams.a1 = 0.0608;
	leafAreaParams.a2 = -0.1293;
	leafAreaParams.b0 = 0.00144;
	leafAreaParams.b1 = 0.0025;
	leafAreaParams.b2 = -0.11;

	// Leaf appearance parameters
	scienceAPI.read("leaf_app_rate1", "", false, appearanceRate1);
	scienceAPI.read("leaf_app_rate2", "", false, appearanceRate2);
	scienceAPI.read("leaf_no_rate_change", "", false, noRateChange);

	scienceAPI.read("tillerSdIntercept", "", false, tillerSdIntercept); // Eqn 13
	scienceAPI.read("tillerSdSlope", "", false, tillerSdSlope); // Eqn 13
	scienceAPI.read("tillerSlaBound", "", false, tillerSlaBound); // Reeves Eqn 
	scienceAPI.read("maxLAIForTillerAddition", "", false, maxLAIForTillerAddition); // Reeves Eqn 






	double temp = -99.0;
	scienceAPI.get("tillerSdIntercept", "", true, temp);
	if (temp > -99)tillerSdIntercept = temp;

	temp = -99.0;
	scienceAPI.get("tillerSdSlope", "", true, temp);
	if (temp > -99)tillerSdSlope = temp;
}
//------------------------------------------------------------------------------------------------
//----------- update Leaf state variables at the end of the day
//------------------------------------------------------------------------------------------------
void LeafCulms::updateVars(void)
{
	leafAppearance.clear();
	culmArea.clear();
	culmLAI.clear();
	tillers = 0.0;
	for (int i = 0; i < (int)Culms.size(); ++i)
	{
		//Culms[i]->updateVars();
		tillers += Culms[i]->getProportion();
		leafAppearance.push_back(Culms[i]->getCurrentLeafNo());
		culmArea.push_back(Culms[i]->getLeafArea());
		culmLAI.push_back(Culms[i]->getTotalLAI());
	}
	tillers--;
	Culms[0]->getCurrentLeafNo();
	Leaf::updateVars();
}

//--------------------------------------------------------------------------------------------------
// Register variables for other modules
//--------------------------------------------------------------------------------------------------
void LeafCulms::doRegistrations(void)
{
	//"expose" only works if the object has been created at the Init1 event
	Leaf::doRegistrations();
	scienceAPI.expose("tillersAdded", "()", "Number of Tillers Added", false, tillersAdded);
	scienceAPI.expose("Supply", "()", "Tiller Supply factor", false, supply);
	scienceAPI.expose("Demand", "()", "Tiller Supply factor", false, demand);
	scienceAPI.expose("calculatedTillers", "()", "Number of Tillers Calculated", false, calculatedTillers);
	scienceAPI.expose("tillers", "()", "Number of Tillers Currently", false, tillers);
	scienceAPI.expose("linearLAI", "()", "Linear LAI", false, linearLAI);
	scienceAPI.expose("tillerSdIntercept", "()", "Linear LAI", false, tillerSdIntercept);
	scienceAPI.expose("tillerSdSlope", "()", "Linear LAI", false, tillerSdSlope);
	scienceAPI.expose("laiReductionForSLA", "()", "LAI Reduction due to carbon limitation", false, laiReductionForSLA);
	scienceAPI.expose("totalLaiReductionForSLA", "()", "Accumulated LAI Reduction due to carbon limitation", false, totalLaiReductionForSLA);
	scienceAPI.expose("maxLaiTarget", "()", "Target LAI value for SLA Target", false, maxLaiTarget);
	scienceAPI.expose("tillerLaiToReduce", "()", "Amount of LAI to reduce to hit SLa Target", false, tillerLaiToReduce);

	scienceAPI.exposeFunction("LeafSizesMain", "mm2", "Size of each leaf on the main culm",
		FloatArrayFunction(&LeafCulms::getLeafSizesMain));
	scienceAPI.exposeFunction("LeafSizesTiller2", "mm2", "Size of each leaf on T2",
		FloatArrayFunction(&LeafCulms::getLeafSizesTiller2));
	scienceAPI.exposeFunction("LeafSizesTiller3", "mm2", "Size of each leaf on T3",
		FloatArrayFunction(&LeafCulms::getLeafSizesTiller3));
	scienceAPI.exposeFunction("LeafSizesTiller4", "mm2", "Size of each leaf on T4",
		FloatArrayFunction(&LeafCulms::getLeafSizesTiller4));
	scienceAPI.exposeFunction("LeafSizesTiller5", "mm2", "Size of each leaf on T5",
		FloatArrayFunction(&LeafCulms::getLeafSizesTiller5));
	scienceAPI.exposeFunction("LeafApp", "()", "Number of leaves on each culm",
		FloatArrayFunction(&LeafCulms::LeafApp));
	scienceAPI.exposeFunction("CulmArea", "()", "Leaf Area on each culm",
		FloatArrayFunction(&LeafCulms::CulmArea));
	scienceAPI.exposeFunction("CulmLAI", "()", "Leaf LAI for each culm",
		FloatArrayFunction(&LeafCulms::CulmLAI));
	scienceAPI.exposeFunction("Proportions", "()", "Leaf LAI for each culm",
		FloatArrayFunction(&LeafCulms::Proportions));
}

void LeafCulms::CulmArea(vector<float>& result)
{
	DVecToFVec(result, culmArea);
}

void LeafCulms::CulmLAI(vector<float>& result)
{
	DVecToFVec(result, culmLAI);
}

void LeafCulms::Proportions(vector<float>& result)
{
	vector<double> props;
	for (int i = 0; i < Culms.size(); ++i)
		props.push_back(Culms[i]->getProportion());
	DVecToFVec(result, props);
}

void LeafCulms::LeafApp(vector<float>& result)
{
	DVecToFVec(result, leafAppearance);
}

void LeafCulms::calcPotentialArea(void)
{
	dltPotentialLAI = 0.0;
	dltStressedLAI = 0.0;

	if (stage > emergence)// && stage <= flag)
	{
		double stressEffect = Min(Min(plant->water->getExpansionStress(), plant->nitrogen->getExpansionStress()), plant->phosphorus->getExpansionStress());
		for (int i = 0; i < (int)Culms.size(); ++i)
		{
			dltStressedLAI += Culms[i]->calcPotentialLeafArea(density, stressEffect);
		}
	}
}

void LeafCulms::areaActual(void)
{
	// Leaves cannot grow too thin - SLA has a maximum point
	laiReductionForSLA = 0.0;
	maxLaiTarget = 0.0;
	tillerLaiToReduce = 0.0;
	SLA = 0.0;
	maxSLA = 0.0;

	laiReductionForSLA = calcCarbonLimitation();

	if (noAreaAdjustmentNeeded())
	{
		dltLAI = dltStressedLAI;
		updateCulmLeafAreas();
		return;
	}

	dltLAI = dltStressedLAI - laiReductionForSLA;
	SLA = calcSLA();

	tillerLaiToReduce = calcCeaseTillerSignal();

	bool moreToAdd = (tillersAdded < calculatedTillers) && (linearLAI < maxLAIForTillerAddition);
	double nLeaves = Culms[0]->getCurrentLeafNo();
	//active tiller reduction is rate limited, so only 0.3 tillers can be deactivated - the rest must come from reducing new leaf growth

	if (nLeaves > 7 && !moreToAdd && tillerLaiToReduce > 0.0)
	{
		double maxTillerLoss = 0.3;				/// externalise as parameter
		double accProportion = 0.0;
		double tillerLaiLeftToReduce = tillerLaiToReduce;

		char msg[120];
		sprintf(msg, "  Cease Tiller. \t\t\tTotal area to remove: %.3f\n", tillerLaiToReduce);
		scienceAPI.write(msg);
		for (unsigned i = Culms.size() - 1; i >= 1; i--)
		{
			if (accProportion < maxTillerLoss && tillerLaiLeftToReduce > 0)
			{
				double tillerArea = Culms[i]->getTotalLAI() + Culms[i]->getDltLAI();
				double tillerProportion = Culms[i]->getProportion();
				if (tillerProportion > 0.0 && tillerArea > 0.0)
				{
					sprintf(msg, "\t Tiller No: %d \t\tProportion: %.3f \t\tArea: %.3f \t\tTiller area to remove: %.3f\n", i, tillerProportion, tillerArea, tillerLaiLeftToReduce);
					scienceAPI.write(msg);

					//use the amount of LAI past the target as an indicator of how much of the tiller
					//to remove which will affect tomorrow's growth - up to the maxTillerLoss
					double propn = Max(0.0, Min(maxTillerLoss - accProportion, tillerLaiLeftToReduce / tillerArea));
					accProportion += propn;
					tillerLaiLeftToReduce -= propn * tillerArea;
					double remainingProportion = Max(0.0, Culms[i]->getProportion() - propn);
					Culms[i]->setProportion(remainingProportion); //can't increase the proportion

					sprintf(msg, "\t Remove Proportion: %.3f \t\tAcc proportion: %.3f \t\tLAI reduction: %.3f, \t\tArea left to remove: %.3f\n\n", propn, accProportion, propn * tillerArea, tillerLaiLeftToReduce);
					scienceAPI.write(msg);

					//if leaf is over sla hard limit, remove as much of the new growth from this tiller first rather than proportionally across all
					double amountToRemove = Min(laiReductionForSLA, Culms[i]->getDltLAI());
					Culms[i]->setDltLAI(Culms[i]->getDltLAI() - amountToRemove);
					laiReductionForSLA -= amountToRemove;
				}
			}
		}
	}

	reduceAllTillersProportionately(laiReductionForSLA);
	updateCulmLeafAreas();
	reportAreaDiscrepency();
}

double LeafCulms::calcCeaseTillerSignal()
{
	// calculate sla target that is below the actual SLA - so as the leaves gets thinner it signals to the tillers to cease growing further
	// max SLA (thinnest leaf) possible using Reeves (1960's Kansas) SLA = 429.72 - 18.158 * LeafNo
	double nLeaves = Culms[0]->getCurrentLeafNo();
	maxSLA = 429.72 - 18.158 * (nLeaves + dltLeafNo);
	maxSLA *= ((100 - tillerSlaBound) / 100.0);		// sla bound vary 30 - 40%
	maxSLA = Min(400, maxSLA);
	maxSLA = Max(150, maxSLA);

	//calc how much LAI we need to remove to get back to the SLA target line
	//this value will be limited by the proportion of tiller area in maxTillerLoss 
	//dltStressedLai can be greater than the actual SLA limitation would allow
	//provides a stronger signal
	maxLaiTarget = maxSLA * (dmGreen + dltDmGreen) / 10000;
	return Max(lai + dltStressedLAI - maxLaiTarget, 0);
}

bool LeafCulms::noAreaAdjustmentNeeded()
{
	//don't calculate SLA limitations prior to endjuv
	if (stage < endJuv)
	{
		updateCulmLeafAreas();
		return true;
	}

	//new leaf growth should have stopped at flagleaf - do we allow for some tillers to finish after flag?
	if (stage >= flag)
	{
		updateCulmLeafAreas();
		return true;
	}
	return false;
}

double LeafCulms::calcCarbonLimitation()
{
	laiReductionForSLA = Max(dltStressedLAI - (dltDmGreen * slaMax * smm2sm), 0.0);
	totalLaiReductionForSLA += laiReductionForSLA;
	if (laiReductionForSLA > 0)
	{
		char msg[120];
		scienceAPI.write(" Leaf Area reduced due to carbon limitation: \n");
		sprintf(msg, "\t dltStressedLAI: %.3f \t\tReduce by: %.3f \t\tdltDmGreen: %.3f\n", dltStressedLAI, laiReductionForSLA, dltDmGreen);
		scienceAPI.write(msg);
	}
	return laiReductionForSLA;
}

double LeafCulms::calcSLA()
{
	if (dmGreen + dltDmGreen <= 0.0) return 0.0;

	return (lai + dltLAI) / (dmGreen + dltDmGreen) * 10000;	// (cm^2/g)
}


void LeafCulms::reportAreaDiscrepency()
{
	double totalDltLeaf = 0.0;
	for (unsigned i = 0; i < Culms.size(); i++) totalDltLeaf += Culms[i]->getDltLAI();

	double diffInLeafArea = totalDltLeaf - dltLAI;
	if (abs(diffInLeafArea) > 0.0001)
	{
		char msg[120];
		scienceAPI.write(" Diff in DltStressedLeaf and Culm Leaf Area Values: \n");
		sprintf(msg, "\t dltStressedLAI: %.5f \t\tTotal DltLAI in Culms: %.5f \t\tDiff: %.7f \n", dltStressedLAI, totalDltLeaf, diffInLeafArea);
		scienceAPI.write(msg);
	}

	double totalLAI = 0.0;
	for (unsigned i = 0; i < Culms.size(); i++) totalLAI += Culms[i]->getTotalLAI();

	diffInLeafArea = totalLAI - lai;
	if (abs(diffInLeafArea) > 0.000001)
	{
		char msg[120];
		scienceAPI.write(" Diff in Leaf LAI and Culm Leaf LAI Values: \n");
		sprintf(msg, "\t LAI: %.3f \t\tTotal LAI in Culms: %.3f \t\tDiff: %.7f \n", lai, totalLAI, diffInLeafArea);
		scienceAPI.write(msg);
	}
}

void LeafCulms::reduceAllTillersProportionately(double laiReduction)
{
	if (laiReduction <= 0.0) return;

	double totalDltLeaf = 0.0;
	for (unsigned i = 0; i < Culms.size(); i++) totalDltLeaf += Culms[i]->getDltLAI();

	//reduce new leaf growth proportionally across all culms
	//not reducing the number of tillers at this stage
	if (totalDltLeaf > 0.0)
	{
		for (unsigned i = 0; i < Culms.size(); i++)
		{
			double dLAI = Culms[i]->getDltLAI();
			//adjust culm dltLAI by proportion of total dltLAI
			double culmProportionToRemove = Max(dLAI / totalDltLeaf * laiReduction, 0);
			Culms[i]->setDltLAI(dLAI - culmProportionToRemove);
		}
	}
}
void LeafCulms::updateCulmLeafAreas()
{
	for (unsigned i = 0; i < Culms.size(); i++)
	{
		Culms[i]->setTotalLAI(Culms[i]->getTotalLAI() + Culms[i]->getDltLAI());
	}
}

void LeafCulms::calcLeafNo(void)
{
	// Overriding this function to retain existing code on Leaf class, but changing this one to manage tillers and leaves.
	// First culm is the main one, that also provides timing for remaining tiller appearance.

	// Calculate final leaf number up to initiation and update each culm FLN and canopy params.
	if (Culms.size() > 0 && stage >= emergence)
	{
		// Update canopy parameters for each tiller.
		if (stage <= fi)
		{
			calcFinalLeafNo();
			for (int i = 0; i < (int)Culms.size(); ++i)
			{
				leafAreaParams.aMaxMain = leafAreaParams.aMaxS * finalLeafNo + leafAreaParams.aMaxI;	// Area of the largest leaf on Main.
				leafAreaParams.X0Main = leafAreaParams.aX0S * finalLeafNo + leafAreaParams.aX0I;		// Posn of largest leaf on Main.
				Culms[i]->setCanopyParams(leafAreaParams, finalLeafNo);
			}
		}

		// 
		// Calculate leaf growth for each culm. Tillers take longer than main to get all leaves expanded.
		// Stop calculating at SGF for now.
		int currentLeaf = (int)floor(Culms[0]->getCurrentLeafNo());
		if (stage < startGrainFill)
		{
			dltLeafNo = Culms[0]->calcLeafAppearance(plant->phenology->getDltTT(), appearanceRate1, appearanceRate2, noRateChange);
			for (int i = 1; i < (int)Culms.size(); ++i)
			{
				Culms[i]->calcLeafAppearance(plant->phenology->getDltTT(), appearanceRate1, appearanceRate2, noRateChange);
			}
		}
		// Calculate tiller numbers
		if (currentLeaf > startThermalQuotientLeafNo) calcTillers(currentLeaf);

	}
}
void LeafCulms::calcTillers(int currentLeaf)
{
	// Calculate tiller addition
	int newLeaf = (int)floor(Culms[0]->getCurrentLeafNo());

	// Up to L5 FE store PTQ. At L5 FE calculate tiller number (endThermalQuotientLeafNo).
	// At L5 FE newLeaf = 6 and currentLeaf = 5
	if (newLeaf >= startThermalQuotientLeafNo + 1 && currentLeaf < endThermalQuotientLeafNo + 1)
	{
		//need to calculate the average R/oCd per day during leaf 5 expansion
		radiationValues += plant->today.radn;
		temperatureValues += plant->phenology->getDltTT();

		if (newLeaf == endThermalQuotientLeafNo + 1)	// L5 Fully Expanded
		{
			double PTQ = radiationValues / temperatureValues;
			calcTillerNumber(PTQ);
			AddInitialTillers();
		}
	}

	// Each time a leaf becomes fully expanded starting at 5 see if a tiller should be initiated.
	// When a leaf is fully expanded a tiller can be initiated at the axil 3 leaves less
	// So at L5 FE (newLeaf = 6, currentLeaf = 5) a Tiller might be at axil 2. i.e. a T2 

	// Add any new tillers and then calc each tiller in turn.
	// Add a tiller if 1. There are more tillers to add.
	//					2. linearLAI < maxLAIForTillerAddition
	//					3. A leaf has fully expanded.  (newLeaf >= 6, newLeaf > currentLeaf)
	//					4. there should be a tiller at that node. ( Check tillerOrder)
	tillersAdded = Culms.size() - 1;
	double lLAI = calcLinearLAI();
	if (newLeaf >= 6 && newLeaf > currentLeaf && calculatedTillers > tillersAdded && calcLinearLAI() < maxLAIForTillerAddition)
	{
		// Axil = currentLeaf - 3
		int newNodeNumber = currentLeaf - 3;
		if (std::find(tillerOrder.begin(), tillerOrder.end(), newNodeNumber) != tillerOrder.end())
		{
			double fractionToAdd = min(1.0, calculatedTillers - tillersAdded);
			initiateTiller(newNodeNumber, fractionToAdd, 1);
		}
	}
}

void LeafCulms::calcTillerNumber(double PTQ)
{

	//the final tiller number (Ftn) is calculated after the full appearance of LeafNo 5 - when leaf 6 emerges.
	//Calc Supply = R/oCd * LA5 * Phy5
	double L5Area = Culms[0]->leafSizes[4];
	double L9Area = Culms[0]->leafSizes[8];
	double Phy5 = appearanceRate1;

	supply = PTQ * L5Area * Phy5;
	//Calc Demand = LA9 - LA5
	demand = L9Area - L5Area;
	double sd = supply / demand;
	calculatedTillers = tillerSdIntercept + tillerSdSlope * sd;
	calculatedTillers = Max(calculatedTillers, 0.0);
	//	calculatedTillers = min(calculatedTillers, 5.0);

	char msg[120];
	sprintf(msg, "Calculated Tiller Number = %.3f\n", calculatedTillers);
	scienceAPI.write(msg);
	sprintf(msg, "Calculated Supply = %.3f\n", supply);
	scienceAPI.write(msg);
	sprintf(msg, "Calculated Demand = %.3f\n", demand);
	scienceAPI.write(msg);

}

void LeafCulms::AddInitialTillers(void)
{
	// OLD COMMENTS
	//tiller emergence is more closely aligned with tip apearance, but we don't track tip, so will use ligule appearance
	//could also use Thermal Time calcs if needed
	//Environmental & Genotypic Control of Tillering in Sorghum ppt - Hae Koo Kim
	//T2=L3, T3=L4, T4=L5, T5=L6

	// If 3 or more then add full T2, initiate a T3 then daily increase tiller number at leaf_initiation_rate. T3, T4�

	//logic to add new tillers depends on which tiller, which is defined by calculatedTillers
	//2 tillers = T3 + T4
	//3 tillers = T2 + T3 + T4
	//4 tillers = T2 + T3 + T4 + T5
	//more than that is too many tillers - but will assume existing pattern for 3 and 4
	//5 tillers = T2 + T3 + T4 + T5 + T6


	// NEW
	// Lafarge et al. (2002) reported a common hierarchy of tiller emergence of T3>T4>T2>T1>T5>T6 across diverse density treatments
	//1 tiller  = T3 
	//2 tillers = T3 + T4
	//3 tillers = T2 + T3 + T4
	//4 tillers = T1 + T2 + T3 + T4
	//5 tillers = T1 + T2 + T3 + T4 + T5
	//6 tillers = T1 + T2 + T3 + T4 + T5 + T6

	// At leaf 5 fully expanded only initialize T1 with 2 leaves if present.

	tillerOrder.clear();
	int nTillers = ceil(calculatedTillers);
	if (nTillers <= 0)return;

	if (nTillers < 3)tillerOrder.push_back(3);
	if (nTillers == 2)tillerOrder.push_back(4);
	if (nTillers == 3)
	{
		tillerOrder.push_back(2);
		tillerOrder.push_back(3);
		tillerOrder.push_back(4);
	}
	if (nTillers > 3)
	{
		for (int i = 1; i <= nTillers; i++)
		{
			tillerOrder.push_back(i);
		}
	}

	if (nTillers > 3)
		initiateTiller(1, 1, 1);
}

void LeafCulms::initiateTiller(int tillerNumber, double fractionToAdd, double initialLeaf)
{
	Culm* newCulm = new Culm(tillerNumber);
	newCulm->setCanopyParams(leafAreaParams, finalLeafNo);
	newCulm->setCurrentLeafNo(initialLeaf);
	newCulm->setProportion(fractionToAdd);

	Culms.push_back(newCulm);
}

/*
void LeafCulms::addTillerProportion(double leafAtAppearance, double fractionToAdd)
{
	//Add a fraction of a tiller every day.
	double currentTillerFraction = Culms.back()->getProportion();

	if (currentTillerFraction + fractionToAdd > 1)
	{
		//update the last tiller to be 1 and add new tiller
		Culms.back()->setProportion(1.0);

		initiateTiller(Culms.back()->getCulmNo() + 1, currentTillerFraction + fractionToAdd - 1.0, 1);
	}
	else
	{
		Culms.back()->setProportion(currentTillerFraction + fractionToAdd);
	}

}*/

/*
void LeafCulms::calcTillerAppearance(int newLeafNo, int currentLeafNo)
{
	//if there are still more tillers to add
	//and the newleaf is greater than 3

	// get number of tillers added so far


	if (calculatedTillers > tillersAdded)
	{
		// calculate linear LAI
		double pltsPerMetre = plant->getPlantDensity() * plant->getRowSpacing() / 1000.0 * plant->getSkipRow();
		linearLAI = pltsPerMetre * tpla / 10000.0;

		if (linearLAI < maxLAIForTillerAddition)
		{
			double fractionToAdd = Min(plant->phenology->getDltTT() / appearanceRate1, calculatedTillers - tillersAdded);
			addTillerProportion(1, fractionToAdd);
			tillersAdded += fractionToAdd;

		}
	}
}*/

double LeafCulms::calcLinearLAI(void)
{
	// calculate linear LAI
	double pltsPerMetre = plant->getPlantDensity() * plant->getRowSpacing() / 1000.0 * plant->getSkipRow();
	linearLAI = pltsPerMetre * tpla / 10000.0;
	return linearLAI;
}

void LeafCulms::getLeafSizesMain(vector<float>& result)
{

	DVecToFVec(result, Culms[0]->leafSizes);
}

void LeafCulms::getLeafSizesTiller2(vector<float>& result)
{
	if (Culms.size() > 1)
	{
		DVecToFVec(result, Culms[1]->leafSizes);
	}
	else
	{
		DVecToFVec(result, vector<double>());
	}
}
void LeafCulms::getLeafSizesTiller3(vector<float>& result)
{
	if (Culms.size() > 2)
	{
		DVecToFVec(result, Culms[2]->leafSizes);
	}
	else
	{
		DVecToFVec(result, vector<double>());
	}
}
void LeafCulms::getLeafSizesTiller4(vector<float>& result)
{
	if (Culms.size() > 3)
	{
		DVecToFVec(result, Culms[3]->leafSizes);
	}
	else
	{
		DVecToFVec(result, vector<double>());
	}
}
void LeafCulms::getLeafSizesTiller5(vector<float>& result)
{
	if (Culms.size() > 4)
	{
		DVecToFVec(result, Culms[4]->leafSizes);
	}
	else
	{
		DVecToFVec(result, vector<double>());
	}
}

//------------------------------------------------------------------------------------------------
//------ LeafCulms_Fixed
//------------------------------------------------------------------------------------------------

//------------------------------------------------------------------------------------------------
//------ LeafCulms_Fixed Constructor
//------------------------------------------------------------------------------------------------
LeafCulms_Fixed::LeafCulms_Fixed(ScienceAPI2& api, Plant* p) : LeafCulms(api, p)
{
}
//------------------------------------------------------------------------------------------------
//------ Destructor
//------------------------------------------------------------------------------------------------
LeafCulms_Fixed::~LeafCulms_Fixed()
{
}

void LeafCulms_Fixed::readParams(void)
{
	Leaf::readParams();
	//for (int i = 0; i < (int)Culms.size(); ++i)
	//	Culms[i]->getParams();

	//scienceAPI.read("aMaxVert", "", false, aMaxVert); // Eqn 13
	//scienceAPI.read("aTillerVert", "", false, aTillerVert); // Eqn 13
}

void LeafCulms_Fixed::calcLeafNo(void)
{
	//overriding this function to retain existing code on Leaf class, but changing this one to manage tillers and leaves
	//first culm is the main one, that also provides timing for remaining tiller appearance

	// calculate final leaf number up to initiation - would finalLeafNo need a different calc for each culm?

	if (Culms.size() > 0 && stage >= emergence)
	{
		//calc finalLeafNo for first culm (main), and then calc it's dltLeafNo
		//add any new tillers and then calc each tiller in turn
		if (stage <= fi)
		{
			//Culms[0]->calcFinalLeafNo();
			//Culms[0]->setCulmNo(0);
			finalLeafNo = Culms[0]->getFinalLeafNo();
			//Culms[0]->calculateLeafSizes();

		}
		double currentLeafNo = Culms[0]->getCurrentLeafNo();
		//double dltLeafNoMainCulm = Culms[0]->calcLeafAppearance();
		//dltLeafNo = dltLeafNoMainCulm; //updates nLeaves
		double newLeafNo = Culms[0]->getCurrentLeafNo();

		calcTillerAppearance((int)floor(newLeafNo), (int)floor(currentLeafNo));

		for (int i = 1; i < (int)Culms.size(); ++i)
		{
			if (stage <= fi)
			{
				//Culms[i]->calcFinalLeafNo();
				//Culms[i]->calculateLeafSizes();
			}
			//Culms[i]->calcLeafAppearance();
		}
	}
}
void LeafCulms_Fixed::calcTillerAppearance(int newLeafNo, int currentLeafNo)
{
	//if there are still more tillers to add
	//and the newleaf is greater than 3
	if (plant->getFtn() > tillersAdded)
	{
		//tiller emergence is more closely aligned with tip apearance, but we don't track tip, so will use ligule appearance
		//could also use Thermal Time calcs if needed
		//Environmental & Genotypic Control of Tillering in Sorghum ppt - Hae Koo Kim
		//T2=L3, T3=L4, T4=L5, T5=L6

		//logic to add new tillers depends on which tiller, which is defined by FTN (fertileTillerNo)
		//this should be provided at sowing  //what if fertileTillers == 1?
		//2 tillers = T3 + T4
		//3 tillers = T2 + T3 + T4
		//4 tillers = T2 + T3 + T4 + T5
		//more than that is too many tillers - but will assume existing pattern for 3 and 4
		//5 tillers = T2 + T3 + T4 + T5 + T6

		bool newLeaf = newLeafNo > currentLeafNo;
		bool newTiller = newLeaf && newLeafNo >= 3; //is it a new leaf, and it is leaf 3 or more
		if (newTiller)
		{
			//tiller 2 emergences with leaf 3, and then adds 1 each time
			//not sure what I'm supposed to do with tiller 1
			//if there are only 2 tillers, then t2 is not present - T3 & T4 are
			//if there is a fraction - between 2 and 3, 
			//this can be interpreted as a proportion of plants that have 2 and a proportion that have 3. 
			//to keep it simple, the fraction will be applied to the 2nd tiller
			double leafAppearance = Culms.size() + 2; //first culm added will equal 3
			double fraction = 1.0;
			if (plant->getFtn() > 2 && plant->getFtn() < 3 && leafAppearance < 4)
			{
				fraction = fmod(plant->getFtn(), 1);
				//tillersAdded += fraction;
			}
			else
			{
				if (plant->getFtn() - tillersAdded < 1)
					fraction = plant->getFtn() - tillersAdded;
				//tillersAdded += 1;
			}

			initiateTiller(Culms.size(), fraction, leafAppearance);
			tillersAdded += fraction;
			////a new tiller is created with each new leaf, up the number of fertileTillers
			//Culm* newCulm = new Culm(scienceAPI, plant, leafAppearance);
			//newCulm->readParams();
			//newCulm->setCurrentLeafNo(leafAppearance-1);
			//verticalAdjustment = aMaxVert;
			//verticalAdjustment += (Culms.size() - 1) * 0.05;
			//newCulm->setVertLeafAdj(verticalAdjustment);
			//newCulm->setProportion(fraction);
			//newCulm->calcFinalLeafNo();
			//newCulm->calcLeafAppearance();
			//newCulm->calculateLeafSizes();
			//Culms.push_back(newCulm);

			//bell curve distribution is adjusted horizontally by moving the curve to the left.
			//This will cause the first leaf to have the same value as the nth leaf on the main culm.
			//T3&T4 were defined during dicussion at initial tillering meeting 27/06/12
			//all others are an assumption
			//T2 = 3 Leaves
			//T3 = 4 Leaves
			//T4 = 5 leaves
			//T5 = 6 leaves
			//T6 = 7 leaves
		}
	}
}
void LeafCulms_Fixed::initiateTiller(int tillerNumber, double fractionToAdd, double initialLeaf)
{
	double leafNoAtAppearance = 1.0;							// DEBUG  parameter?
	double nTillersPresent = Culms.size() - 1;
	/*
	Culm* newCulm = new Culm(scienceAPI, plant, leafNoAtAppearance);
	newCulm->getParams();
	newCulm->setCulmNo(tillerNumber);
	newCulm->setCurrentLeafNo(initialLeaf);
	verticalAdjustment = aMaxVert + aTillerVert * (nTillersPresent - 1);
	newCulm->setVertLeafAdj(verticalAdjustment);
	newCulm->setProportion(fractionToAdd);
	newCulm->calcFinalLeafNo();

	//	newCulm->calcLeafAppearance();
	newCulm->calculateLeafSizes();
	Culms.push_back(newCulm);
	*/
}

void LeafCulms_Fixed::areaActual(void)
{
	dltLAI = dltStressedLAI;
	if (stage >= endJuv && stage < flag)
	{
		laiReductionForSLA = calcCarbonLimitation();
		dltLAI = dltStressedLAI - laiReductionForSLA;
		SLA = calcSLA();
	}
}
