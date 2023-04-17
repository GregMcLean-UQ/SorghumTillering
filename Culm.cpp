//------------------------------------------------------------------------------------------------
#include <stdio.h>
#include "Plant.h"
#include "LeafCulms.h"
#include "Culm.h"


using namespace Sorghum;

//------------------------------------------------------------------------------------------------
//------ Culm Constructor
//------------------------------------------------------------------------------------------------
Culm::Culm()
{
	
	initialize();
	
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
	// leaf number
	finalLeafNo = 0.0;
	dltLeafNo = 0.0;
	largestLeafPlateau = 0; //value less than 1 means ignore it
	currentLeafNo = 1.0;//noEmergence - is in Leaf
	//finalLeafCorrection = 0;
	//vertAdjValue = 0.0;
	proportion = 1.0;
	totalLAI = 0.0;
	dltStressedLAI = 0.0;
	dltLAI = 0.0;
	culmNo = 0;
	leafArea = 0;
	totalArea = 0;

	//readParams();
}


//------------------------------------------------------------------------------------------------
//-----------  read leaf parameters
//------------------------------------------------------------------------------------------------
void Culm::initT(LeafAreaParams leafAreaParams, int posn)
{

	leafNoCorrection = leafAreaParams.leafNoCorrection;       //corrects for other growing leaves
	largestLeafPlateau = leafAreaParams.largestLeafPlateau;

	culmNo = posn;

}
void Culm::readParams(void)
{

 }
void Culm::setCanopyParams(LeafAreaParams lAP, double finalLeafNumber)
{
	finalLeafNo = finalLeafNumber - culmNo;
	// Canopy parameters
	// aMax and x0 for tillers have values relative to those of the main stem.
	// Relative size factor = 0, 0.23, 0.13, 0.13, 0.13, 0.39 ...
	double relLeafSize = 0;
	if (culmNo == 0)relLeafSize = 0.0;
	else if (culmNo == 1)relLeafSize = 0.23;
	else if (culmNo < 5)relLeafSize = 0.13;
	else relLeafSize = 0.39;

	x0 = lAP.X0Main - (culmNo == 0 ? 0 : culmNo + 1);
	aMax = lAP.aMaxMain * (1 - relLeafSize);

	a = lAP.a0 + (lAP.a1 / (1 + lAP.a2 * finalLeafNo));	// Bell shaped curve characteristics.
	b = lAP.b0 + (lAP.b1 / (1 + lAP.b2 * finalLeafNo));
	leafSizes.clear();
	for (int leafNo = 1; leafNo < ceil(finalLeafNo) + 1; leafNo++)
		leafSizes.push_back(aMax * exp(a * pow((leafNo - x0), 2) + b * pow((leafNo - x0), 3)) * 100);
}

double Culm::calcLeafAppearance(double dltTT, double appearanceRate1, double  appearanceRate2, double noRateChange)
{
	// Calculate increase in leaf expandion for today.
	dltLeafNo = 0.0;
	double remainingLeaves = finalLeafNo - currentLeafNo;
	if (remainingLeaves > 0.0)
	{
		// Peter's 2 stage version used here, modified to apply to last few leaves before flag
		// i.e. c_leaf_no_rate_change is leaf number from the top down (e.g. 4)
		double leafAppRate = (remainingLeaves <= noRateChange) ? appearanceRate2 : appearanceRate1;

		// if leaves are still growing, the cumulative number of phyllochrons or fully expanded
		// leaves is calculated from thermal time for the day.
		dltLeafNo = bound(divide(dltTT, leafAppRate), 0.0, remainingLeaves);
		currentLeafNo += dltLeafNo;
	}
	return dltLeafNo;
}


double Culm::calcPotentialLeafArea(void)
{
	// Today's total culm area is calculated. 
	double leafNoEffective = Min(currentLeafNo + leafNoCorrection, finalLeafNo);
	double leafAreaNow = leafArea;
	leafArea = 0.0;
	for (int i = 0; i < ceil(leafNoEffective); i++)
	{
		double fraction = min(1.0, leafNoEffective - i);
		leafArea += leafSizes[i] * fraction;
	}
	return max(leafArea - leafAreaNow,0.0);
}


/*	double leafsize = calcIndividualLeafSize(leafNoEffective);
	//leafArea = getAreaOfCurrentLeaf(leafNoEffective);		HACK
	//leafArea *= proportion; //proportion is 1 unless this tiller is a fraction ie: Fertile Tiller Number is 2.2, then 1 tiller is 0.2
	totalArea += leafsize * dltLeafNo;
	leafArea = leafsize * smm2sm * density * dltLeafNo; // in dltLai
	dltLAI = leafArea * proportion;
	return dltLAI;
	//totalLAI += leafArea;
	//return (leafArea * proportion);
}*/

double Culm::getAreaOfCurrentLeaf(double leafNo)
{
	// interpolate leaf sizes to get area of this leaf
	// check upper
	if (leafNo > leafSizes.size())
		return leafSizes.back();
	else
	{
		int leafIndx = (int)floor(leafNo) - 1;
		double leafPart = leafNo - floor(leafNo);
		double size = leafSizes[leafIndx] + (leafSizes[leafIndx + 1] - leafSizes[leafIndx]) * leafPart;
		return size;
	}
}


double Culm::calcIndividualLeafSize(double leafNo)
{
	// use finalLeafNo to calculate the size of the individual leafs
	// Eqn 5 from Improved methods for predicting individual leaf area and leaf senescence in maize
	// (Zea mays) C.J. Birch, G.L. Hammer and K.G. Ricket. Aust. J Agric. Res., 1998, 49, 249-62
	//


	//double leafPlateauStart = 24;
	//adding new code to handle varieties that grow very high number of leaves
	/*

	if (largestLeafPlateau > 1)
	{
		if (finalLeafNo > largestLeafPlateau)
		{
			largestLeafPos = aX0I + aX0S * largestLeafPlateau;

			if (leafNo > largestLeafPos)
			{
				double tailCount = largestLeafPlateau - largestLeafPos;
				if (leafNo < finalLeafNo - tailCount)
				{
					leafNo = largestLeafPos;
				}
				else
				{
					leafNo = largestLeafPlateau - (finalLeafNo - leafNo);
				}
			}
		}
	}*/

	//Relationship for calculating maximum individual leaf area from Total Leaf No
	//Source: Modelling genotypic and environmental control of leaf area dynamics in grain sorghum. II. Individual leaf level 
	//Carberry, Muchow, Hammer,1992
	//written as Y = Y0*exp(a*pow(X-X0,2)+b*(pow(X-X0,3))) 
	//pg314 -Leaf area production model

	//Largest Leaf calculation
	//originally from "Improved methods for predicting individual leaf area and leaf senescence in maize" - Birch, Hammer, Rickert 1998
	//double aMaxB = 4.629148, aMaxC = 6.6261562; 
	//double aMax = aMaxA * (1 - exp(-aMaxB * (finalLeafNo - aMaxC)));  // maximum individual leaf area
	//Calculation then changed to use the relationship as described in the Carberry paper in Table 2
	//The actual intercept and slope will be determined by the cultivar, and read from the config file (sorghum.xml)
	//aMaxS = 19.5; //not 100% sure what this number should be - tried a range and this provided the best fit forthe test data

	//double leafSize = largestLeafSize * exp(a * pow((leafNo - largestLeafPos), 2) + b * pow((leafNo - largestLeafPos), 3)) * 100;
	return 0;// leafSize;
}

