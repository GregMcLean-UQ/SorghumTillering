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
	class Culm //: public PlantComponent
	{
		// private Methods -------------------------------------------------------
	private:

		//double leafNoAtAppearance;

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
		double dltStressedLAI; //adjusted growth for stress for the current day


		//double lastLeafNumber;
		//double finalLeafCorrection;
		//double vertAdjValue;

		
		int culmNo;

		// plant
		double density;
	public:
		vector<double> leafSizes;

		// public Methods -------------------------------------------------------
	public:
		Culm(ScienceAPI2& api, Plant* p, double leafNoAtApp);
		~Culm();

		
		///virtual void initialize(void);
		 //void initT(LeafAreaParams leafAreaParams, int posn);
		//virtual void updateVars(void);
		//virtual void doRegistrations(void);

		 void initialize(void) ;
		 void readParams(void);
		 //void updateVars(void) ;
		void initT(LeafAreaParams leafAreaParams, int posn);
		 void doRegistrations(void);

		void setCanopyParams(LeafAreaParams lAP, double finalLeafNumber);

		double getCurrentLeafNo(void);
		void setCurrentLeafNo(const double& val);
		//double getlastLeafNumber() { return lastLeafNumber; }
		double gettotalArea() { return totalArea; }

		//void setLastLeafNo(double _leafNo) { lastLeafNumber = _leafNo; }
		double getFinalLeafNo(void);
		double calcLeafAppearance(double dltTT, double appearanceRate1, double  appearanceRate2, double noRateChange);
		double calcPotentialLeafArea(void);
		double calcIndividualLeafSize(double leafNo);
		void setProportion(double proportion);
		double getLeafArea();
		double getProportion() { return proportion; }

		void setCulmNo(int _culmNo) { culmNo = _culmNo; }
		int getCulmNo() { return culmNo; }
		double getAreaOfCurrentLeaf(double leaves);

		double getTotalLAI() { return totalLAI; }
		void setTotalLAI(double val) { totalLAI = val; }
		double getDltLAI() { return dltLAI; }
		double setDltLAI(double val) { dltLAI = val; }
		double getStressedDltLAI() { return dltStressedLAI; }
		void setStressedDltLAI(double val) { dltStressedLAI = val; }

	};


}
#endif