#ifndef ADAPPARAM_H
#define ADAPPARAM_H
#include <vector>
#include <fstream>
using namespace std;

class AdapParam{
public:
  AdapParam() {}
  ~AdapParam() {}

  static void initPriesAdapParam();
  static void initRandomAdapParam();
  static void setPara2PSO(double X[]);

  enum OptCate{
	NO_OPT=0,
	PSO=1,
	DOWNHILL=2,
  };
  // 优化方法
  enum OptType{
    STDPSO=0,
	YSPSO=1,
	SELPSO=2,
	QUAPSO=3,
  };

  enum ErrType{
    NO_ADAP_ERR=0,
    FLOWDIR=-1,
    PORDER=-2,
    NORDER=-3,
    NEGDIAM=-4,
    LINEQU=-5,
    NO_CONV=-6,
    ERR_HEMO=-7,
  };

  // 需优化的参数 for wall adaptation
  // 参数[0]:参数值 [1]:上限 [2]:下限
  // PO2 reference value
  static double PO2Ref[3];
  // Tau reference value
  static double tauRef[3];
  // Sigma reference value
  static double sigmaRef[3];
  // w reference value
  static double wRef[3];
  // k coefficients
  static double kc[3], kmd[3], kmg[3], ksd[3], ksg[3], kwtau[3], kwsigma[3];
  // J0
  static double J0[3];
  // Lref
  static double LRef[3];

  // For no wall adaptation
  static double kc2[3], kp2[3], km2[3], ks2[3];  // For NoWallAdap
  static double J02[3], tauRef2[3], QRef2[3], LRef2[3];

  // 其它自适应计算相关的参数
  static double dt;
  static double QRef;

  // Optimization type
  static int optCate;
  static int optType;
  static double ErrorD, ErrorV, ErrorQ;
  static int errFlag;
  static int AdapErrCnt;

  // PSO parameters
  static int nPar;

  // Calculation order
  static vector<int> PosCalcOrder;
  static vector<int> NegCalcOrder;

  // Output stream
  static ofstream adapLogFile;
  static ofstream adapHemoFile;
  static ofstream adapGBestFile;
};

#endif
