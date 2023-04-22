//------------------------------------------------------------------------------------------------
//---------------------------------------------------------------------------

#ifndef CulmH
#define CulmH

#include <vector>
#include <string>
#include "Utilities.h"
#include <iostream>
#include <fstream>
#include <sstream>
#include "LeafCulms.h"
using namespace std;

namespace Sorghum {

	//------------------------------------------------------------------------------------------------
	//------ Culm
	//------------------------------------------------------------------------------------------------
	struct LeafAreaParams;
	class Culm 
	{
		// private Methods -------------------------------------------------------
	private:

		int culmNo;

		// Leaf Area
		double aMax;	//Area of the largest leaf.
		double x0;		// position of the largest leaf.
		double a, b;	// Bell shaped curve characteristics.

		double leafNoCorrection;       //corrects for other growing leaves
		double largestLeafPlateau;
		double finalLeafNo;
		double dltLeafNo;
		double currentLeafNo;
		double proportion;
		double leafArea;//changes each day
		double totalLAI;		// accumulated lai for this culm
		double totalArea;	// total area of the culm
		double dltLAI; //growth for the current day

	public:
		vector<double> leafSizes;

		// public Methods -------------------------------------------------------
	public:
		Culm(int posn);
		~Culm();

		void initialize(void);

		void setCanopyParams(LeafAreaParams lAP, double finalLeafNumber);

		double gettotalArea() { return totalArea; }

		double calcLeafAppearance(double dltTT, double appearanceRate1, double  appearanceRate2, double noRateChange);
		double calcPotentialLeafArea(double density, double stressEffect);
		double calcIndividualLeafSize(double leafNo);
		double getProportion() { return proportion; }

		int getCulmNo() { return culmNo; }

		double Culm::culmArea(double nLeaves);

		double getTotalLAI() { return totalLAI; }
		void setTotalLAI(double val) { totalLAI = val; }
		double getDltLAI() { return dltLAI; }
		void setDltLAI(double val) { dltLAI = val; }
		void setProportion(double val) { proportion = val; }
		double getLeafArea() { return leafArea; }

		double getFinalLeafNo(void) { return finalLeafNo; }
		double getCurrentLeafNo(void) { return currentLeafNo; }
		void setCurrentLeafNo(const double& val) { currentLeafNo = val; }

	};


}
#endif