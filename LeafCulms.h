//------------------------------------------------------------------------------------------------
//---------------------------------------------------------------------------

#ifndef LeafCulmsH
#define LeafCulmsH

#include <vector>
#include <string>
#include "Utilities.h"
#include <iostream>
#include <fstream>
#include <sstream>
#include "Culm.h"

using namespace std;

namespace Sorghum {
	//------------------------------------------------------------------------------------------------
	struct LeafAreaParams
	{
		//Birch, Hammer bell shaped curve parameters
		double aX0I;                    // Eqn 14 calc x0 - position of largest leaf
		double aX0S;                    
		double X0Main;
		double aMaxI;
		double aMaxS;
		double aMaxMain;
		double a0, a1, a2;
		double b0, b1, b2;
		double leafNoCorrection;
		double largestLeafPlateau;
	};
	class LeafCulms : public Leaf
	{
		// private Methods -------------------------------------------------------
	private:
		LeafAreaParams leafAreaParams;
	protected:
		vector<Culm*> Culms;
		vector<int> tillerOrder;
//		double verticalAdjustment;
		double calculatedTillers;
		double supply;
		double demand;

		double radiationValues;
		double temperatureValues;
		double avgRadiation;
		double thermalTimeCount;
		double maxLAIForTillerAddition;
		int startThermalQuotientLeafNo;
		int endThermalQuotientLeafNo;

		double tillerSdIntercept;
		double tillerSdSlope;
		double tillerSlaBound;
		double linearLAI;
		double laiReductionForSLA;
		double totalLaiReductionForSLA;
		double maxLaiTarget;
		double tillerLaiToReduce;

		// leaf appearance
		double noSeed;
		double noEmergence;
		double initRate;
		//double initialTPLA;
		//double tplaInflectionRatio,tplaProductionCoef;
		
		double appearanceRate1;
		double appearanceRate2;
		double noRateChange;
		double minLeafNo, maxLeafNo;
		// public Methods -------------------------------------------------------
	public:
		LeafCulms(ScienceAPI2&, Plant* p);
		virtual ~LeafCulms();

		virtual void  initialize(void);
		virtual void  readParams(void);
		virtual void  updateVars(void);
		virtual void  doRegistrations(void);

		virtual void calcLeafNo(void);
		virtual void calcPotentialArea(void);
		virtual void areaActual(void);

		void calcTillers(int currentLeaf);
		//virtual void calcTillerAppearance(int newLeafNo, int currentLeafNo);
		void calcTillerNumber(double PTQ);
		void AddInitialTillers(void);
		virtual void initiateTiller(int tillerNumber, double fractionToAdd, double initialLeaves);
		double calcLinearLAI(void);
		void LeafCulms::calculateTillerCessation(void);


		//void addTillerProportion(double leafAtAppearance, double fractionToAdd);
		virtual double calcCeaseTillerSignal();
		virtual double calcCarbonLimitation();
		virtual double calcSLA();
		virtual void reportAreaDiscrepency();
		virtual void reduceAllTillersProportionately(double laiReduction);
		virtual void updateCulmLeafAreas();

		void getLeafSizesMain(vector<float>& result);
		void getLeafSizesTiller1(vector<float>& result);
		void getLeafSizesTiller2(vector<float>& result);
		void getLeafSizesTiller3(vector<float>& result);
		void getLeafSizesTiller4(vector<float>& result);
		void getLeafSizesTiller5(vector<float>& result);
		void LeafApp(vector<float>& result);
		void CulmArea(vector<float>& result);
		void CulmLAI(vector<float>& result);
		void Proportions(vector<float>& result);

		vector<double> leafAppearance;
		vector<double> culmArea;
		vector<double> culmLAI;

		virtual bool allFullyExpanded(void);

	};

	class LeafCulms_Fixed : public LeafCulms
	{
		// private Methods -------------------------------------------------------
	private:

		// public Methods -------------------------------------------------------
	public:
		LeafCulms_Fixed(ScienceAPI2&, Plant* p);
		virtual ~LeafCulms_Fixed();

		virtual void readParams(void);
		virtual void calcLeafNo(void);
		virtual void calcTillerAppearance(int newLeafNo, int currentLeafNo);
		virtual void areaActual(void);
		virtual void initiateTiller(int tillerNumber, double fractionToAdd, double initialLeaf);
	};
	//------------------------------------------------------------------------------------------------
}
#endif
