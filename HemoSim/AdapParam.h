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
	static void initYJLAdapParam();
	static void initAdapParamBounder();
	static void setPara2PSO(double X[]);

	enum OptCate{
		NO_OPT=0,
		PSO=1,
		DOWNHILL=2,
		PSA=3, /*参数灵敏度分析*/
		RST=4, /*结果分析*/
	};
	// 优化方法
	enum OptType{
		STDPSO=0, /*例子中粒子数目20，学习因子都取2，惯性权重0.7，迭代步数取10000*/
		CMPPSO=1, /*学习因子相加大于4，典型：c1=2.8,c2=1.3 粒子数目30，学习因子都为1，迭代步数10000*/
		STDQPSO=2, /*例子中粒子数目20，线性压缩扩张系数，迭代步数取5000*/
		SELQPSO=3, /*例子中粒子数目20，随机压缩扩张系数，带选择算子，迭代步数取5000*/
		SECPSO=4, /*本例中学习因子不能都取2，参考例子给出：粒子数40，c1=c2=1, 惯性权重0.7， 迭代步数10000*/
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
	static int optMethod;
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
