//------------------------------------------------------------------------------------------------
#include <stdio.h>
#include "Plant.h"
#include "LeafCulms.h"
#include "Culm.h"


using namespace Sorghum;

//------------------------------------------------------------------------------------------------
//------ Culm Constructor
//------------------------------------------------------------------------------------------------
Culm::Culm(int posn)
{
	initialize();
	culmNo = posn;
}
//------------------------------------------------------------------------------------------------
//------ Destructor
//------------------------------------------------------------------------------------------------
Culm::~Culm()
{
}
//------------------------------------------------------------------------------------------------
//------- Initialize variables
//------------------------------------------------------------------------------------------------
void Culm::initialize(void)
{
	// Leaf number
	finalLeafNo = 0.0;
	dltLeafNo = 0.0;
	largestLeafPlateau = 0; // Value less than 1 means ignore it.
	currentLeafNo = 1.0;	// noEmergence.
	// Leaf Area
	totalLAI = 0.0;
	dltLAI = 0.0;
	leafArea = 0.0;
	// Tiller proportion
	proportion = 1.0;
}

//------------------------------------------------------------------------------------------------
//-----------  read leaf parameters
//------------------------------------------------------------------------------------------------
void Culm::setCanopyParams(LeafAreaParams lAP, double finalLeafNumber)
{
	// Initialise the Leaf Area parameters.
	finalLeafNo = finalLeafNumber - culmNo;
	// Canopy parameters
	// aMax and x0 for tillers have values relative to those of the main stem.
	// Relative size factor = 0, 0.23, 0.13, 0.13, 0.13, 0.39 ...
	double relLeafSize = 0;
	if (culmNo == 0)relLeafSize = 0.0;
	else if (culmNo == 1)relLeafSize = 0.23;
	else if (culmNo < 5)relLeafSize = 0.13;
	else relLeafSize = 0.39;

	// Calculate the position of the largest leaf (x0) and the size of the largest leaf (aMax).
	x0 = lAP.X0Main - (culmNo == 0 ? 0 : culmNo + 1);
	aMax = lAP.aMaxMain * (1 - relLeafSize);

	// Bell shaped curve characteristics.
	a = lAP.a0 + (lAP.a1 / (1 + lAP.a2 * finalLeafNo));	
	b = lAP.b0 + (lAP.b1 / (1 + lAP.b2 * finalLeafNo));

	// Calculate individual leaf size.
	leafSizes.clear();
	for (int leafNo = 1; leafNo < ceil(finalLeafNo) + 1; leafNo++)
		leafSizes.push_back(aMax * exp(a * pow((leafNo - x0), 2) + b * pow((leafNo - x0), 3)) * 100);
	leafNoCorrection = lAP.leafNoCorrection;       //corrects for other growing leaves
	largestLeafPlateau = lAP.largestLeafPlateau;
}

//------------------------------------------------------------------------------------------------
//-----------  Calculate leaf appearance
//------------------------------------------------------------------------------------------------

double Culm::calcLeafAppearance(double dltTT, double appearanceRate1, double  appearanceRate2, double noRateChange)
{
	// Calculate increase in leaf expansion for today.
	dltLeafNo = 0.0;
	double remainingLeaves = finalLeafNo - currentLeafNo;
	if (remainingLeaves > 0.0)
	{
		// Peter's 2 stage version used here, modified to apply to last few leaves before flag
		// i.e. c_leaf_no_rate_change is leaf number from the top down (e.g. 4)
		double leafAppRate = (remainingLeaves <= noRateChange) ? appearanceRate2 : appearanceRate1;

		// If leaves are still growing, the cumulative number of phyllochrons or fully expanded
		// leaves is calculated from thermal time for the day.
		dltLeafNo = bound(divide(dltTT, leafAppRate), 0.0, remainingLeaves);
		currentLeafNo += dltLeafNo;
	}
	return dltLeafNo;
}

//------------------------------------------------------------------------------------------------
//-----------  Calculate Leaf size and area
//------------------------------------------------------------------------------------------------
double Culm::calcPotentialLeafArea(double density, double stressEffect)
{
	// Today's total culm area is calculated in LAI so  * smm2sm * density reduced by stress
	double leafNoEffective = Min(currentLeafNo + leafNoCorrection, finalLeafNo);
	// Get the potential increase in leaf area due the the leaf expansion for today. 
	double leafAreaNow = calcCulmArea(leafNoEffective - dltLeafNo);
	leafArea = calcCulmArea(leafNoEffective);
	double dltLeafArea = max(leafArea - leafAreaNow, 0.0);
	dltLAI = dltLeafArea * smm2sm * density * stressEffect * proportion;
	return dltLAI;
}

double Culm::calcCulmArea(double nLeaves)
{
	// Sum the area of each leaf plus the fraction of the last.
	double area = 0;
	for (int i = 0; i < ceil(nLeaves); i++)
	{
		double fraction = min(1.0, nLeaves - i);
		area += leafSizes[i] * fraction;
	}
	return area;
}

