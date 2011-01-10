#include "MantidKernel/MagneticIon.h"
#include <sstream>
#include <stdexcept>

namespace Mantid
{
namespace PhysicalConstants
{
// Empty constructor 
MagneticIon::MagneticIon():symbol(""),charge(0)        
{
  std::vector <double> zero(8,0.);
  j0=zero;
  j2=zero;
  j4=zero;
  j6=zero;
}
// constructor from symbol, charge, and 4 vectors of 8 doubles each
MagneticIon::MagneticIon(const std::string symbol, const uint16_t charge, const double j0i[8],
  const double j2i[8],const double j4i[8],const double j6i[8]):symbol(symbol),charge(charge)
{
  
  std::vector<double> temp0(j0i,j0i+8);
  j0=temp0;
  std::vector<double> temp2(j2i,j2i+8);
  j2=temp2;
  std::vector<double> temp4(j4i,j4i+8);
  j4=temp4;
  std::vector<double> temp6(j6i,j6i+8);
  j6=temp6;     
}

// copy constructor
MagneticIon::MagneticIon(const MagneticIon& other): symbol(other.symbol), charge(other.charge), j0(other.j0), j2(other.j2), j4(other.j4), j6(other.j6)
{
}

// initialize all magnetic ions, and store names into ion_map
int initializeMap()
{
  if((int)ion_map.size()<10)
  {
    static const MagneticIon Sc0((std::string )"Sc", (uint16_t)0, (double[]){0.2512,90.03,0.329,39.402,0.4235,14.322,-0.0043,0.2029}, // <j0>
        (double[]){10.8172,54.327,4.7353,14.847,0.6071,4.218,0.0011,0.1212}, // <j2>
        (double[]){1.342,10.2,0.3837,3.079,0.0468,0.118,-0.0328,0.1343}, // <j4>
        (double[]){0,0,0,0,0,0,0,0}); // <j6>
    ion_map["Sc0"]=Sc0;
    static const MagneticIon Sc1((std::string )"Sc", (uint16_t)1, (double[]){0.4889,51.16,0.5203,14.076,-0.0286,0.179,0.0185,0.1217}, // <j0>
        (double[]){8.5021,34.285,3.2116,10.994,0.4244,3.605,0.0009,0.1037}, // <j2>
        (double[]){7.1167,15.487,-6.6671,18.269,0.49,2.992,0.0047,0.1624}, // <j4>
        (double[]){0,0,0,0,0,0,0,0}); // <j6>
    ion_map["Sc1"]=Sc1;
    static const MagneticIon Sc2((std::string )"Sc", (uint16_t)2, (double[]){0.5048,31.403,0.5186,10.99,-0.0241,1.183,0,0.0578}, // <j0>
        (double[]){4.3683,28.654,3.7231,10.823,0.6074,3.668,0.0014,0.0681}, // <j2>
        (double[]){-1.6684,15.648,1.7742,9.062,0.4075,2.412,0.0042,0.1105}, // <j4>
        (double[]){0,0,0,0,0,0,0,0}); // <j6>
    ion_map["Sc2"]=Sc2;
    static const MagneticIon Ti0((std::string )"Ti", (uint16_t)0, (double[]){0.4657,33.59,0.549,9.879,-0.0291,0.323,0.0123,0.1088}, // <j0>
        (double[]){4.3583,36.056,3.823,11.133,0.6855,3.469,0.002,0.0967}, // <j2>
        (double[]){-2.1515,11.271,2.5149,8.859,0.3555,2.149,0.0045,0.1244}, // <j4>
        (double[]){0,0,0,0,0,0,0,0}); // <j6>
    ion_map["Ti0"]=Ti0;
    static const MagneticIon Ti1((std::string )"Ti", (uint16_t)1, (double[]){0.5093,36.703,0.5032,10.371,-0.0263,0.311,0.0116,0.1125}, // <j0>
        (double[]){6.1567,27.275,2.6833,8.983,0.407,3.052,0.0011,0.0902}, // <j2>
        (double[]){-1.0383,16.19,1.4699,8.924,0.3631,2.283,0.0044,0.127}, // <j4>
        (double[]){0,0,0,0,0,0,0,0}); // <j6>
    ion_map["Ti1"]=Ti1;
    static const MagneticIon Ti2((std::string )"Ti", (uint16_t)2, (double[]){0.5091,24.976,0.5162,8.757,-0.0281,0.916,0.0015,0.0589}, // <j0>
        (double[]){4.3107,18.348,2.096,6.797,0.2984,2.548,0.0007,0.064}, // <j2>
        (double[]){-1.3242,15.31,1.2042,7.899,0.3976,2.156,0.0051,0.082}, // <j4>
        (double[]){0,0,0,0,0,0,0,0}); // <j6>
    ion_map["Ti2"]=Ti2;
    static const MagneticIon Ti3((std::string )"Ti", (uint16_t)3, (double[]){0.3571,22.841,0.6688,8.931,-0.0354,0.483,0.0099,0.0575}, // <j0>
        (double[]){3.3717,14.444,1.8258,5.713,0.247,2.265,0.0005,0.0491}, // <j2>
        (double[]){-1.1117,14.635,0.7689,6.927,0.4385,2.089,0.006,0.0572}, // <j4>
        (double[]){0,0,0,0,0,0,0,0}); // <j6>
    ion_map["Ti3"]=Ti3;
    static const MagneticIon V0((std::string )"V", (uint16_t)0, (double[]){0.4086,28.811,0.6077,8.544,-0.0295,0.277,0.0123,0.097}, // <j0>
        (double[]){3.76,21.831,2.4026,7.546,0.4464,2.663,0.0017,0.0556}, // <j2>
        (double[]){-0.9633,15.273,0.9274,7.732,0.3891,2.053,0.0063,0.084}, // <j4>
        (double[]){0,0,0,0,0,0,0,0}); // <j6>
    ion_map["V0"]=V0;
    static const MagneticIon V1((std::string )"V", (uint16_t)1, (double[]){0.4444,32.648,0.5683,9.097,-0.2285,0.022,0.215,0.1111}, // <j0>
        (double[]){4.7474,23.323,2.3609,7.808,0.4105,2.706,0.0014,0.08}, // <j2>
        (double[]){-0.9606,15.545,1.1278,8.118,0.3653,2.097,0.0056,0.1027}, // <j4>
        (double[]){0,0,0,0,0,0,0,0}); // <j6>
    ion_map["V1"]=V1;
    static const MagneticIon V2((std::string )"V", (uint16_t)2, (double[]){0.4085,23.853,0.6091,8.246,-0.1676,0.041,0.1496,0.0593}, // <j0>
        (double[]){3.4386,16.53,1.9638,6.141,0.2997,2.267,0.0009,0.0565}, // <j2>
        (double[]){-1.1729,14.973,0.9092,7.613,0.4105,2.039,0.0067,0.0719}, // <j4>
        (double[]){0,0,0,0,0,0,0,0}); // <j6>
    ion_map["V2"]=V2;
    static const MagneticIon V3((std::string )"V", (uint16_t)3, (double[]){0.3598,19.336,0.6632,7.617,-0.3064,0.03,0.2835,0.0515}, // <j0>
        (double[]){2.3005,14.682,2.0364,6.13,0.4099,2.382,0.0014,0.0252}, // <j2>
        (double[]){-0.9417,14.205,0.5284,6.607,0.4411,1.967,0.0076,0.0569}, // <j4>
        (double[]){0,0,0,0,0,0,0,0}); // <j6>
    ion_map["V3"]=V3;
    static const MagneticIon V4((std::string )"V", (uint16_t)4, (double[]){0.3106,16.816,0.7198,7.049,-0.0521,0.302,0.0221,0.0433}, // <j0>
        (double[]){1.8377,12.267,1.8247,5.458,0.3979,2.248,0.0012,0.0399}, // <j2>
        (double[]){-0.7654,13.097,0.3071,5.674,0.4476,1.871,0.0081,0.0518}, // <j4>
        (double[]){0,0,0,0,0,0,0,0}); // <j6>
    ion_map["V4"]=V4;
    static const MagneticIon Cr0((std::string )"Cr", (uint16_t)0, (double[]){0.1135,45.199,0.3481,19.493,0.5477,7.354,-0.0092,0.1975}, // <j0>
        (double[]){3.4085,20.127,2.1006,6.802,0.4266,2.394,0.0019,0.0662}, // <j2>
        (double[]){-0.667,19.613,0.5342,6.478,0.3641,1.905,0.0073,0.0628}, // <j4>
        (double[]){0,0,0,0,0,0,0,0}); // <j6>
    ion_map["Cr0"]=Cr0;
    static const MagneticIon Cr1((std::string )"Cr", (uint16_t)1, (double[]){-0.0977,0.047,0.4544,26.005,0.5579,7.489,0.0831,0.1114}, // <j0>
        (double[]){3.7768,20.346,2.1028,6.893,0.401,2.411,0.0017,0.0686}, // <j2>
        (double[]){-0.8309,18.043,0.7252,7.531,0.3828,2.003,0.0073,0.0781}, // <j4>
        (double[]){0,0,0,0,0,0,0,0}); // <j6>
    ion_map["Cr1"]=Cr1;
    static const MagneticIon Cr2((std::string )"Cr", (uint16_t)2, (double[]){1.2024,-0.005,0.4158,20.548,0.6032,6.956,-1.2218,0.0572}, // <j0>
        (double[]){2.6422,16.06,1.9198,6.253,0.4446,2.372,0.002,0.048}, // <j2>
        (double[]){-0.893,15.664,0.559,7.033,0.4093,1.924,0.0081,0.0631}, // <j4>
        (double[]){0,0,0,0,0,0,0,0}); // <j6>
    ion_map["Cr2"]=Cr2;
    static const MagneticIon Cr3((std::string )"Cr", (uint16_t)3, (double[]){-0.3094,0.027,0.368,17.035,0.6559,6.524,0.2856,0.0436}, // <j0>
        (double[]){1.6262,15.066,2.0618,6.284,0.5281,2.368,0.0023,0.0263}, // <j2>
        (double[]){-0.7327,14.073,0.3268,5.674,0.4114,1.81,0.0085,0.0505}, // <j4>
        (double[]){0,0,0,0,0,0,0,0}); // <j6>
    ion_map["Cr3"]=Cr3;
    static const MagneticIon Cr4((std::string )"Cr", (uint16_t)4, (double[]){-0.232,0.043,0.3101,14.952,0.7182,6.173,0.2042,0.0419}, // <j0>
        (double[]){1.0293,13.95,1.9933,6.059,0.5974,2.346,0.0027,0.0366}, // <j2>
        (double[]){-0.6748,12.946,0.1805,6.753,0.4526,1.8,0.0098,0.0644}, // <j4>
        (double[]){0,0,0,0,0,0,0,0}); // <j6>
    ion_map["Cr4"]=Cr4;
    static const MagneticIon Mn0((std::string )"Mn", (uint16_t)0, (double[]){0.2438,24.963,0.1472,15.673,0.6189,6.54,-0.0105,0.1748}, // <j0>
        (double[]){2.6681,16.06,1.7561,5.64,0.3675,2.049,0.0017,0.0595}, // <j2>
        (double[]){-0.5452,15.471,0.4406,4.902,0.2884,1.543,0.0059,0.0488}, // <j4>
        (double[]){0,0,0,0,0,0,0,0}); // <j6>
    ion_map["Mn0"]=Mn0;
    static const MagneticIon Mn1((std::string )"Mn", (uint16_t)1, (double[]){-0.0138,0.421,0.4231,24.668,0.5905,6.655,-0.001,0.1242}, // <j0>
        (double[]){3.2953,18.695,1.8792,6.24,0.3927,2.201,0.0022,0.0659}, // <j2>
        (double[]){-0.7947,17.867,0.6078,7.704,0.3798,1.905,0.0087,0.0737}, // <j4>
        (double[]){0,0,0,0,0,0,0,0}); // <j6>
    ion_map["Mn1"]=Mn1;
    static const MagneticIon Mn2((std::string )"Mn", (uint16_t)2, (double[]){0.422,17.684,0.5948,6.005,0.0043,-0.609,-0.0219,0.0589}, // <j0>
        (double[]){2.0515,15.556,1.8841,6.063,0.4787,2.232,0.0027,0.0306}, // <j2>
        (double[]){-0.7416,15.255,0.3831,6.469,0.3935,1.8,0.0093,0.0577}, // <j4>
        (double[]){0,0,0,0,0,0,0,0}); // <j6>
    ion_map["Mn2"]=Mn2;
    static const MagneticIon Mn3((std::string )"Mn", (uint16_t)3, (double[]){0.4198,14.283,0.6054,5.469,0.9241,-0.009,-0.9498,0.0392}, // <j0>
        (double[]){1.2427,14.997,1.9567,6.118,0.5732,2.258,0.0031,0.0336}, // <j2>
        (double[]){-0.6603,13.607,0.2322,6.218,0.4104,1.74,0.0101,0.0579}, // <j4>
        (double[]){0,0,0,0,0,0,0,0}); // <j6>
    ion_map["Mn3"]=Mn3;
    static const MagneticIon Mn4((std::string )"Mn", (uint16_t)4, (double[]){0.376,12.566,0.6602,5.133,-0.0372,0.563,0.0011,0.0393}, // <j0>
        (double[]){0.7879,13.886,1.8717,5.743,0.5981,2.182,0.0034,0.0434}, // <j2>
        (double[]){-0.5127,13.461,0.0313,7.763,0.4282,1.701,0.0113,0.0693}, // <j4>
        (double[]){0,0,0,0,0,0,0,0}); // <j6>
    ion_map["Mn4"]=Mn4;
    static const MagneticIon Fe0((std::string )"Fe", (uint16_t)0, (double[]){0.0706,35.008,0.3589,15.358,0.5819,5.561,-0.0114,0.1398}, // <j0>
        (double[]){1.9405,18.473,1.9566,6.323,0.5166,2.161,0.0036,0.0394}, // <j2>
        (double[]){-0.5029,19.677,0.2999,3.776,0.2576,1.424,0.0071,0.0292}, // <j4>
        (double[]){0,0,0,0,0,0,0,0}); // <j6>
    ion_map["Fe0"]=Fe0;
    static const MagneticIon Fe1((std::string )"Fe", (uint16_t)1, (double[]){0.1251,34.963,0.3629,15.514,0.5223,5.591,-0.0105,0.1301}, // <j0>
        (double[]){2.629,18.66,1.8704,6.331,0.469,2.163,0.0031,0.0491}, // <j2>
        (double[]){-0.5109,19.25,0.3896,4.891,0.281,1.526,0.0069,0.0375}, // <j4>
        (double[]){0,0,0,0,0,0,0,0}); // <j6>
    ion_map["Fe1"]=Fe1;
    static const MagneticIon Fe2((std::string )"Fe", (uint16_t)2, (double[]){0.0263,34.96,0.3668,15.943,0.6188,5.594,-0.0119,0.1437}, // <j0>
        (double[]){1.649,16.559,1.9064,6.133,0.5206,2.137,0.0035,0.0335}, // <j2>
        (double[]){-0.5401,17.227,0.2865,3.742,0.2658,1.424,0.0076,0.0278}, // <j4>
        (double[]){0,0,0,0,0,0,0,0}); // <j6>
    ion_map["Fe2"]=Fe2;
    static const MagneticIon Fe3((std::string )"Fe", (uint16_t)3, (double[]){0.3972,13.244,0.6295,4.903,-0.0314,0.35,0.0044,0.0441}, // <j0>
        (double[]){1.3602,11.998,1.5188,5.003,0.4705,1.991,0.0038,0.0374}, // <j2>
        (double[]){-0.5507,11.493,0.2153,4.906,0.3468,1.523,0.0095,0.0314}, // <j4>
        (double[]){0,0,0,0,0,0,0,0}); // <j6>
    ion_map["Fe3"]=Fe3;
    static const MagneticIon Fe4((std::string )"Fe", (uint16_t)4, (double[]){0.3782,11.38,0.6556,4.592,-0.0346,0.483,0.0005,0.0362}, // <j0>
        (double[]){1.5582,8.275,1.1863,3.279,0.1366,1.107,-0.0022,0.0327}, // <j2>
        (double[]){-0.5352,9.507,0.1783,5.175,0.3584,1.469,0.0097,0.036}, // <j4>
        (double[]){0,0,0,0,0,0,0,0}); // <j6>
    ion_map["Fe4"]=Fe4;
    static const MagneticIon Co0((std::string )"Co", (uint16_t)0, (double[]){0.4139,16.162,0.6013,4.78,-0.1518,0.021,0.1345,0.1033}, // <j0>
        (double[]){1.9678,14.17,1.4911,4.948,0.3844,1.797,0.0027,0.0452}, // <j2>
        (double[]){-0.4221,14.195,0.29,3.979,0.2469,1.286,0.0063,0.04}, // <j4>
        (double[]){0,0,0,0,0,0,0,0}); // <j6>
    ion_map["Co0"]=Co0;
    static const MagneticIon Co1((std::string )"Co", (uint16_t)1, (double[]){0.099,33.125,0.3645,15.177,0.547,5.008,-0.0109,0.0983}, // <j0>
        (double[]){2.4097,16.161,1.578,5.46,0.4095,1.914,0.0031,0.0581}, // <j2>
        (double[]){-0.4115,14.561,0.358,4.717,0.2644,1.418,0.0074,0.0541}, // <j4>
        (double[]){0,0,0,0,0,0,0,0}); // <j6>
    ion_map["Co1"]=Co1;
    static const MagneticIon Co2((std::string )"Co", (uint16_t)2, (double[]){0.4332,14.355,0.5857,4.608,-0.0382,0.134,0.0179,0.0711}, // <j0>
        (double[]){1.9049,11.644,1.3159,4.357,0.3146,1.645,0.0017,0.0459}, // <j2>
        (double[]){0.4759,14.046,0.2747,3.731,0.2458,1.25,0.0057,0.0282}, // <j4>
        (double[]){0,0,0,0,0,0,0,0}); // <j6>
    ion_map["Co2"]=Co2;
    static const MagneticIon Co3((std::string )"Co", (uint16_t)3, (double[]){0.3902,12.508,0.6324,4.457,-0.15,0.034,0.1272,0.0515}, // <j0>
        (double[]){1.7058,8.859,1.1409,3.309,0.1474,1.09,-0.0025,0.0462}, // <j2>
        (double[]){-0.4466,13.391,0.1419,3.011,0.2773,1.335,0.0093,0.0341}, // <j4>
        (double[]){0,0,0,0,0,0,0,0}); // <j6>
    ion_map["Co3"]=Co3;
    static const MagneticIon Co4((std::string )"Co", (uint16_t)4, (double[]){0.3515,10.778,0.6778,4.234,-0.0389,0.241,0.0098,0.039}, // <j0>
        (double[]){1.311,8.025,1.1551,3.179,0.1608,1.13,-0.0011,0.0374}, // <j2>
        (double[]){-0.4091,13.194,-0.0194,3.417,0.3534,1.421,0.0112,0.0622}, // <j4>
        (double[]){0,0,0,0,0,0,0,0}); // <j6>
    ion_map["Co4"]=Co4;
    static const MagneticIon Ni0((std::string )"Ni", (uint16_t)0, (double[]){-0.0172,35.739,0.3174,14.269,0.7136,4.566,-0.0143,0.1072}, // <j0>
        (double[]){1.0302,12.252,1.4669,4.745,0.4521,1.744,0.0036,0.0338}, // <j2>
        (double[]){-0.4428,14.485,0.087,3.234,0.2932,1.331,0.0096,0.0554}, // <j4>
        (double[]){0,0,0,0,0,0,0,0}); // <j6>
    ion_map["Ni0"]=Ni0;
    static const MagneticIon Ni1((std::string )"Ni", (uint16_t)1, (double[]){0.0705,35.856,0.3984,13.804,0.5427,4.397,-0.0118,0.0738}, // <j0>
        (double[]){2.104,14.866,1.4302,5.071,0.4031,1.778,0.0034,0.0561}, // <j2>
        (double[]){-0.3836,13.425,0.3116,4.462,0.2471,1.309,0.0079,0.0515}, // <j4>
        (double[]){0,0,0,0,0,0,0,0}); // <j6>
    ion_map["Ni1"]=Ni1;
    static const MagneticIon Ni2((std::string )"Ni", (uint16_t)2, (double[]){0.0163,35.883,0.3916,13.223,0.6052,4.339,-0.0133,0.0817}, // <j0>
        (double[]){1.708,11.016,1.2147,4.103,0.315,1.533,0.0018,0.0446}, // <j2>
        (double[]){-0.3803,10.403,0.2838,3.378,0.2108,1.104,0.005,0.0474}, // <j4>
        (double[]){0,0,0,0,0,0,0,0}); // <j6>
    ion_map["Ni2"]=Ni2;
    static const MagneticIon Ni3((std::string )"Ni", (uint16_t)3, (double[]){0.0012,35,0.3468,11.987,0.6667,4.252,-0.0148,0.0883}, // <j0>
        (double[]){1.4683,8.671,0.1794,1.106,1.1068,3.257,-0.0023,0.0373}, // <j2>
        (double[]){-0.4014,9.046,0.2314,3.075,0.2192,1.084,0.006,0.0323}, // <j4>
        (double[]){0,0,0,0,0,0,0,0}); // <j6>
    ion_map["Ni3"]=Ni3;
    static const MagneticIon Ni4((std::string )"Ni", (uint16_t)4, (double[]){-0.009,35.861,0.2776,11.79,0.7474,4.201,-0.0163,0.0966}, // <j0>
        (double[]){1.1612,7.7,1.0027,3.263,0.2719,1.378,0.0025,0.0326}, // <j2>
        (double[]){-0.3509,8.157,0.222,2.106,0.1567,0.925,0.0065,0.0352}, // <j4>
        (double[]){0,0,0,0,0,0,0,0}); // <j6>
    ion_map["Ni4"]=Ni4;
    static const MagneticIon Cu0((std::string )"Cu", (uint16_t)0, (double[]){0.0909,34.984,0.4088,11.443,0.5128,3.825,-0.0124,0.0513}, // <j0>
        (double[]){1.9182,14.49,1.3329,4.73,0.3842,1.639,0.0035,0.0617}, // <j2>
        (double[]){-0.3204,15.132,0.2335,4.021,0.2312,1.196,0.0068,0.0457}, // <j4>
        (double[]){0,0,0,0,0,0,0,0}); // <j6>
    ion_map["Cu0"]=Cu0;
    static const MagneticIon Cu1((std::string )"Cu", (uint16_t)1, (double[]){0.0749,34.966,0.4147,11.764,0.5238,3.85,-0.0127,0.0591}, // <j0>
        (double[]){1.8814,13.433,1.2809,4.545,0.3646,1.602,0.0033,0.059}, // <j2>
        (double[]){-0.3572,15.125,0.2336,3.966,0.2315,1.197,0.007,0.0397}, // <j4>
        (double[]){0,0,0,0,0,0,0,0}); // <j6>
    ion_map["Cu1"]=Cu1;
    static const MagneticIon Cu2((std::string )"Cu", (uint16_t)2, (double[]){0.0232,34.969,0.4023,11.564,0.5882,3.843,-0.0137,0.0532}, // <j0>
        (double[]){1.5189,10.478,1.1512,3.813,0.2918,1.398,0.0017,0.0429}, // <j2>
        (double[]){-0.3914,14.74,0.1275,3.384,0.2548,1.255,0.0103,0.0394}, // <j4>
        (double[]){0,0,0,0,0,0,0,0}); // <j6>
    ion_map["Cu2"]=Cu2;
    static const MagneticIon Cu3((std::string )"Cu", (uint16_t)3, (double[]){0.0031,34.907,0.3582,10.914,0.6531,3.828,-0.0147,0.0665}, // <j0>
        (double[]){1.2797,8.45,1.0315,3.28,0.2401,1.25,0.0015,0.0389}, // <j2>
        (double[]){-0.3671,14.082,-0.0078,3.315,0.3154,1.377,0.0132,0.0534}, // <j4>
        (double[]){0,0,0,0,0,0,0,0}); // <j6>
    ion_map["Cu3"]=Cu3;
    static const MagneticIon Cu4((std::string )"Cu", (uint16_t)4, (double[]){-0.0132,30.682,0.2801,11.163,0.749,3.817,-0.0165,0.0767}, // <j0>
        (double[]){0.9568,7.448,0.9099,3.396,0.3729,1.494,0.0049,0.033}, // <j2>
        (double[]){-0.2915,14.124,-0.1065,4.201,0.3247,1.352,0.0148,0.0579}, // <j4>
        (double[]){0,0,0,0,0,0,0,0}); // <j6>
    ion_map["Cu4"]=Cu4;
    static const MagneticIon Y0((std::string )"Y", (uint16_t)0, (double[]){0.5915,67.608,1.5123,17.9,-1.113,14.136,0.008,0.3272}, // <j0>
        (double[]){14.4084,44.658,5.1045,14.904,-0.0535,3.319,0.0028,0.1093}, // <j2>
        (double[]){-8.0767,32.201,7.9197,25.156,1.4067,6.827,-0.0001,0.1031}, // <j4>
        (double[]){0,0,0,0,0,0,0,0}); // <j6>
    ion_map["Y0"]=Y0;
    static const MagneticIon Zr0((std::string )"Zr", (uint16_t)0, (double[]){0.4106,59.996,1.0543,18.648,-0.4751,10.54,0.0106,0.3667}, // <j0>
        (double[]){10.1378,35.337,4.7734,12.545,-0.0489,2.672,0.0036,0.0912}, // <j2>
        (double[]){-5.2697,32.868,4.193,24.183,1.5202,6.048,-0.0002,0.0855}, // <j4>
        (double[]){0,0,0,0,0,0,0,0}); // <j6>
    ion_map["Zr0"]=Zr0;
    static const MagneticIon Zr1((std::string )"Zr", (uint16_t)1, (double[]){0.4532,59.595,0.7834,21.436,-0.2451,9.036,0.0098,0.3639}, // <j0>
        (double[]){11.8722,34.92,4.0502,12.127,-0.0632,2.828,0.0034,0.0737}, // <j2>
        (double[]){-5.6384,33.607,4.6729,22.338,1.3258,5.924,-0.0003,0.0674}, // <j4>
        (double[]){0,0,0,0,0,0,0,0}); // <j6>
    ion_map["Zr1"]=Zr1;
    static const MagneticIon Nb0((std::string )"Nb", (uint16_t)0, (double[]){0.3946,49.23,1.3197,14.822,-0.7269,9.616,0.0129,0.3659}, // <j0>
        (double[]){7.4796,33.179,5.0884,11.571,-0.0281,1.564,0.0047,0.0944}, // <j2>
        (double[]){-3.1377,25.595,2.3411,16.569,1.2304,4.99,-0.0005,0.0615}, // <j4>
        (double[]){0,0,0,0,0,0,0,0}); // <j6>
    ion_map["Nb0"]=Nb0;
    static const MagneticIon Nb1((std::string )"Nb", (uint16_t)1, (double[]){0.4572,49.918,1.0274,15.726,-0.4962,9.157,0.0118,0.3403}, // <j0>
        (double[]){8.7735,33.285,4.6556,11.605,-0.0268,1.539,0.0044,0.0855}, // <j2>
        (double[]){-3.3598,25.82,2.8297,16.427,1.1203,4.982,-0.0005,0.0724}, // <j4>
        (double[]){0,0,0,0,0,0,0,0}); // <j6>
    ion_map["Nb1"]=Nb1;
    static const MagneticIon Mo0((std::string )"Mo", (uint16_t)0, (double[]){0.1806,49.057,1.2306,14.786,-0.4268,6.987,0.0171,0.4135}, // <j0>
        (double[]){5.118,23.422,4.1809,9.208,-0.0505,1.743,0.0053,0.0655}, // <j2>
        (double[]){-2.886,20.572,1.813,14.628,1.1899,4.264,-0.0008,0.041}, // <j4>
        (double[]){0,0,0,0,0,0,0,0}); // <j6>
    ion_map["Mo0"]=Mo0;
    static const MagneticIon Mo1((std::string )"Mo", (uint16_t)1, (double[]){0.35,48.035,1.0305,15.06,-0.3929,7.479,0.0139,0.351}, // <j0>
        (double[]){7.2367,28.128,4.0705,9.923,-0.0317,1.455,0.0049,0.0798}, // <j2>
        (double[]){-3.2618,25.486,2.3596,16.462,1.1164,4.491,-0.0007,0.0592}, // <j4>
        (double[]){0,0,0,0,0,0,0,0}); // <j6>
    ion_map["Mo1"]=Mo1;
    static const MagneticIon Tc0((std::string )"Tc", (uint16_t)0, (double[]){0.1298,49.661,1.1656,14.131,-0.3134,5.513,0.0195,0.3869}, // <j0>
        (double[]){4.2441,21.397,3.9439,8.375,-0.0371,1.187,0.0066,0.0645}, // <j2>
        (double[]){-2.7975,20.159,1.652,16.261,1.1726,3.943,-0.0008,0.0657}, // <j4>
        (double[]){0,0,0,0,0,0,0,0}); // <j6>
    ion_map["Tc0"]=Tc0;
    static const MagneticIon Tc1((std::string )"Tc", (uint16_t)1, (double[]){0.2674,48.957,0.9569,15.141,-0.2387,5.458,0.016,0.3412}, // <j0>
        (double[]){6.4056,24.824,3.54,8.611,-0.0366,1.485,0.0044,0.0806}, // <j2>
        (double[]){-2.047,19.683,1.6306,11.592,0.8698,3.769,-0.001,0.0723}, // <j4>
        (double[]){0,0,0,0,0,0,0,0}); // <j6>
    ion_map["Tc1"]=Tc1;
    static const MagneticIon Ru0((std::string )"Ru", (uint16_t)0, (double[]){0.1069,49.424,1.1912,12.742,-0.3176,4.912,0.0213,0.3597}, // <j0>
        (double[]){3.7445,18.613,3.4749,7.42,-0.0363,1.007,0.0073,0.0533}, // <j2>
        (double[]){-1.5042,17.949,0.6027,9.961,0.97,3.393,-0.001,0.0338}, // <j4>
        (double[]){0,0,0,0,0,0,0,0}); // <j6>
    ion_map["Ru0"]=Ru0;
    static const MagneticIon Ru1((std::string )"Ru", (uint16_t)1, (double[]){0.441,33.309,1.4775,9.553,-0.9361,6.722,0.0176,0.2608}, // <j0>
        (double[]){5.2826,23.683,3.5813,8.152,-0.0257,0.426,0.0131,0.083}, // <j2>
        (double[]){1.6278,18.506,1.1828,10.189,0.8138,3.418,-0.0009,0.0673}, // <j4>
        (double[]){0,0,0,0,0,0,0,0}); // <j6>
    ion_map["Ru1"]=Ru1;
    static const MagneticIon Rh0((std::string )"Rh", (uint16_t)0, (double[]){0.0976,49.882,1.1601,11.831,-0.2789,4.127,0.0234,0.3263}, // <j0>
        (double[]){3.3651,17.344,3.2121,6.804,-0.035,0.503,0.0146,0.0545}, // <j2>
        (double[]){-1.3492,17.577,0.4527,10.507,0.9285,3.155,-0.0009,0.0483}, // <j4>
        (double[]){0,0,0,0,0,0,0,0}); // <j6>
    ion_map["Rh0"]=Rh0;
    static const MagneticIon Rh1((std::string )"Rh", (uint16_t)1, (double[]){0.3342,29.756,1.2209,9.438,-0.5755,5.332,0.021,0.2574}, // <j0>
        (double[]){4.026,18.95,3.1663,7,-0.0296,0.486,0.0127,0.0629}, // <j2>
        (double[]){-1.4673,17.957,0.7381,9.944,0.8485,3.126,-0.0012,0.0487}, // <j4>
        (double[]){0,0,0,0,0,0,0,0}); // <j6>
    ion_map["Rh1"]=Rh1;
    static const MagneticIon Pd0((std::string )"Pd", (uint16_t)0, (double[]){0.2003,29.363,1.1446,9.599,-0.3689,4.042,0.0251,0.2453}, // <j0>
        (double[]){3.3105,14.726,2.6332,5.862,-0.0437,1.13,0.0053,0.0492}, // <j2>
        (double[]){-1.1955,17.628,0.3183,11.309,0.8696,2.909,-0.0006,0.0555}, // <j4>
        (double[]){0,0,0,0,0,0,0,0}); // <j6>
    ion_map["Pd0"]=Pd0;
    static const MagneticIon Pd1((std::string )"Pd", (uint16_t)1, (double[]){0.5033,24.504,1.9982,6.908,-1.524,5.513,0.0213,0.1909}, // <j0>
        (double[]){4.2749,17.9,2.7021,6.354,-0.0258,0.7,0.0071,0.0768}, // <j2>
        (double[]){-1.4098,17.765,0.7927,9.999,0.771,2.93,-0.0006,0.053}, // <j4>
        (double[]){0,0,0,0,0,0,0,0}); // <j6>
    ion_map["Pd1"]=Pd1;
    static const MagneticIon Ce2((std::string )"Ce", (uint16_t)2, (double[]){0.2953,17.685,0.2923,6.733,0.4313,5.383,-0.0194,0.0845}, // <j0>
        (double[]){0.9809,18.063,1.8413,7.769,0.9905,2.845,0.012,0.0448}, // <j2>
        (double[]){-0.6468,10.533,0.4052,5.624,0.3412,1.535,0.008,0.0522}, // <j4>
        (double[]){-0.1212,7.994,-0.0639,4.024,0.1519,1.096,0.0078,0.0388}); // <j6>
    ion_map["Ce2"]=Ce2;
    static const MagneticIon Nd2((std::string )"Nd", (uint16_t)2, (double[]){0.1645,25.045,0.2522,11.978,0.6012,4.946,-0.018,0.0668}, // <j0>
        (double[]){1.453,18.34,1.6196,7.285,0.8752,2.622,0.0126,0.0461}, // <j2>
        (double[]){-0.5416,12.204,0.3571,6.169,0.3154,1.485,0.0098,0.0519}, // <j4>
        (double[]){-0.16,8.009,0.0272,4.028,0.1104,1.068,0.0139,0.0363}); // <j6>
    ion_map["Nd2"]=Nd2;
    static const MagneticIon Nd3((std::string )"Nd", (uint16_t)3, (double[]){0.054,25.029,0.3101,12.102,0.6575,4.722,-0.0216,0.0478}, // <j0>
        (double[]){0.6751,18.342,1.6272,7.26,0.9644,2.602,0.015,0.045}, // <j2>
        (double[]){-0.4053,14.014,0.0329,7.005,0.3759,1.707,0.0209,0.0372}, // <j4>
        (double[]){0.0416,8.014,-0.1261,4.04,0.14,1.087,0.0102,0.0367}); // <j6>
    ion_map["Nd3"]=Nd3;
    static const MagneticIon Sm2((std::string )"Sm", (uint16_t)2, (double[]){0.0909,25.203,0.3037,11.856,0.625,4.237,-0.02,0.0408}, // <j0>
        (double[]){1.036,18.425,1.4769,7.032,0.881,2.437,0.0152,0.0345}, // <j2>
        (double[]){-0.415,14.057,0.1368,7.032,0.3272,1.582,0.0192,0.0319}, // <j4>
        (double[]){0.1428,6.041,0.0723,2.033,0.055,0.513,0.0081,0.045}); // <j6>
    ion_map["Sm2"]=Sm2;
    static const MagneticIon Sm3((std::string )"Sm", (uint16_t)3, (double[]){0.0288,25.207,0.2973,11.831,0.6954,4.212,-0.0213,0.051}, // <j0>
        (double[]){0.4707,18.43,1.4261,7.034,0.9574,2.439,0.0182,0.051}, // <j2>
        (double[]){-0.4288,10.052,0.1782,5.019,0.2833,1.236,0.0088,0.0328}, // <j4>
        (double[]){-0.0944,6.03,-0.0498,2.074,0.1372,0.645,-0.0132,0.0387}); // <j6>
    ion_map["Sm3"]=Sm3;
    static const MagneticIon Eu2((std::string )"Eu", (uint16_t)2, (double[]){0.0755,25.296,0.3001,11.599,0.6438,4.025,-0.0196,0.0488}, // <j0>
        (double[]){0.897,18.443,1.3769,7.005,0.906,2.421,0.019,0.0511}, // <j2>
        (double[]){-0.4145,10.193,0.2447,5.164,0.2661,1.205,0.0065,0.0516}, // <j4>
        (double[]){-0.1252,6.049,0.0507,2.085,0.0572,0.646,0.0132,0.0403}); // <j6>
    ion_map["Eu2"]=Eu2;
    static const MagneticIon Eu3((std::string )"Eu", (uint16_t)3, (double[]){0.0204,25.308,0.301,11.474,0.7005,3.942,-0.022,0.0356}, // <j0>
        (double[]){0.3985,18.451,1.3307,6.956,0.9603,2.378,0.0197,0.0447}, // <j2>
        (double[]){-0.4095,10.211,0.1485,5.175,0.272,1.237,0.0131,0.0494}, // <j4>
        (double[]){-0.0817,6.039,-0.0596,2.12,0.1243,0.764,-0.0001,0.0206}); // <j6>
    ion_map["Eu3"]=Eu3;
    static const MagneticIon Gd2((std::string )"Gd", (uint16_t)2, (double[]){0.0636,25.382,0.3033,11.212,0.6528,3.788,-0.0199,0.0486}, // <j0>
        (double[]){0.7756,18.469,1.3124,6.899,0.8956,2.338,0.0199,0.0441}, // <j2>
        (double[]){-0.3824,10.344,0.1955,5.306,0.2622,1.203,0.0097,0.0363}, // <j4>
        (double[]){-0.1351,5.03,0.0828,2.025,0.0315,0.503,0.0187,0.0453}); // <j6>
    ion_map["Gd2"]=Gd2;
    static const MagneticIon Gd3((std::string )"Gd", (uint16_t)3, (double[]){0.0186,25.387,0.2895,11.142,0.7135,3.752,-0.0217,0.0489}, // <j0>
        (double[]){0.3347,18.476,1.2465,6.877,0.9537,2.318,0.0217,0.0484}, // <j2>
        (double[]){-0.3621,10.353,0.1016,5.31,0.2649,1.219,0.0147,0.0494}, // <j4>
        (double[]){-0.0662,6.031,-0.085,2.154,0.1323,0.891,0.0048,0.0371}); // <j6>
    ion_map["Gd3"]=Gd3;
    static const MagneticIon Tb2((std::string )"Tb", (uint16_t)2, (double[]){0.0547,25.509,0.3171,10.591,0.649,3.517,-0.0212,0.0342}, // <j0>
        (double[]){0.6688,18.491,1.2487,6.822,0.8888,2.275,0.0215,0.0439}, // <j2>
        (double[]){-0.3443,10.469,0.1481,5.416,0.2575,1.182,0.0104,0.028}, // <j4>
        (double[]){-0.0758,6.032,-0.054,2.158,0.1199,0.89,0.0051,0.0488}); // <j6>
    ion_map["Tb2"]=Tb2;
    static const MagneticIon Tb3((std::string )"Tb", (uint16_t)3, (double[]){0.0177,25.51,0.2921,10.577,0.7133,3.512,-0.0231,0.0512}, // <j0>
        (double[]){0.2892,18.497,1.1678,6.797,0.9437,2.257,0.0232,0.0458}, // <j2>
        (double[]){-0.3228,10.476,0.0638,5.419,0.2566,1.196,0.0159,0.0439}, // <j4>
        (double[]){-0.0559,6.031,-0.102,2.237,0.1264,1.107,0.0167,0.017}); // <j6>
    ion_map["Tb3"]=Tb3;
    static const MagneticIon Dy2((std::string )"Dy", (uint16_t)2, (double[]){0.1308,18.316,0.3118,7.665,0.5795,3.147,-0.0226,0.0315}, // <j0>
        (double[]){0.5917,18.511,1.1828,6.747,0.8801,2.214,0.0229,0.0439}, // <j2>
        (double[]){-0.3206,12.071,0.0904,8.026,0.2616,1.23,0.0143,0.0767}, // <j4>
        (double[]){-0.0568,6.032,-0.1003,2.24,0.1401,1.106,0.0109,0.0463}); // <j6>
    ion_map["Dy2"]=Dy2;
    static const MagneticIon Dy3((std::string )"Dy", (uint16_t)3, (double[]){0.1157,15.073,0.327,6.799,0.5821,3.02,-0.0249,0.0146}, // <j0>
        (double[]){0.2523,18.517,1.0914,6.736,0.9345,2.208,0.025,0.0476}, // <j2>
        (double[]){-0.2829,9.525,0.0565,4.429,0.2437,1.066,0.0092,0.0181}, // <j4>
        (double[]){-0.0423,6.038,-0.1248,2.244,0.1359,1.2,0.0188,0.035}); // <j6>
    ion_map["Dy3"]=Dy3;
    static const MagneticIon Ho2((std::string )"Ho", (uint16_t)2, (double[]){0.0995,18.176,0.3305,7.856,0.5921,2.98,-0.023,0.124}, // <j0>
        (double[]){0.5094,18.515,1.1234,6.706,0.8727,2.159,0.0242,0.056}, // <j2>
        (double[]){-0.2976,9.719,0.1224,4.635,0.2279,1.005,0.0063,0.0452}, // <j4>
        (double[]){-0.0725,6.045,-0.0318,2.243,0.0738,1.202,0.0252,0.0634}); // <j6>
    ion_map["Ho2"]=Ho2;
    static const MagneticIon Ho3((std::string )"Ho", (uint16_t)3, (double[]){0.0566,18.318,0.3365,7.688,0.6317,2.943,-0.0248,0.0068}, // <j0>
        (double[]){0.2188,18.516,1.024,6.707,0.9251,2.161,0.0268,0.0503}, // <j2>
        (double[]){-0.2717,9.731,0.0474,4.638,0.2292,1.047,0.0124,0.031}, // <j4>
        (double[]){-0.0289,6.05,-0.1545,2.23,0.155,1.26,0.0177,0.0351}); // <j6>
    ion_map["Ho3"]=Ho3;
    static const MagneticIon Er2((std::string )"Er", (uint16_t)2, (double[]){0.1122,18.122,0.3462,6.911,0.5649,2.761,-0.0235,0.0207}, // <j0>
        (double[]){0.4693,18.528,1.0545,6.649,0.8679,2.12,0.0261,0.0413}, // <j2>
        (double[]){-0.2975,9.829,0.1189,4.741,0.2116,1.004,0.0117,0.0524}, // <j4>
        (double[]){0.0648,6.056,-0.0515,2.23,0.0825,1.264,0.025,0.0409}); // <j6>
    ion_map["Er2"]=Er2;
    static const MagneticIon Er3((std::string )"Er", (uint16_t)3, (double[]){0.0586,17.98,0.354,7.096,0.6126,2.748,-0.0251,0.0171}, // <j0>
        (double[]){0.171,18.534,0.9879,6.625,0.9044,2.1,0.0278,0.0489}, // <j2>
        (double[]){-0.2568,9.834,0.0356,4.741,0.2172,1.028,0.0148,0.0434}, // <j4>
        (double[]){-0.011,6.061,-0.1954,2.224,0.1818,1.296,0.0149,0.0455}); // <j6>
    ion_map["Er3"]=Er3;
    static const MagneticIon Tm2((std::string )"Tm", (uint16_t)2, (double[]){0.0983,18.324,0.338,6.918,0.5875,2.662,-0.0241,0.0404}, // <j0>
        (double[]){0.4198,18.542,0.9959,6.6,0.8593,2.082,0.0284,0.0457}, // <j2>
        (double[]){-0.2677,9.888,0.0925,4.784,0.2056,0.99,0.0124,0.0396}, // <j4>
        (double[]){0.0842,4.07,0.0807,0.849,-0.2087,0.039,0.2095,0.036}); // <j6>
    ion_map["Tm2"]=Tm2;
    static const MagneticIon Tm3((std::string )"Tm", (uint16_t)3, (double[]){0.0581,15.092,0.2787,7.801,0.6854,2.793,-0.0224,0.0351}, // <j0>
        (double[]){0.176,18.542,0.9105,6.579,0.897,2.062,0.0294,0.0468}, // <j2>
        (double[]){-0.2292,9.895,0.0124,4.785,0.2108,1.007,0.0151,0.0334}, // <j4>
        (double[]){0.0727,4.073,0.0243,0.689,3.9459,0.002,-3.9076,0.0502}); // <j6>
    ion_map["Tm3"]=Tm3;
    static const MagneticIon Yb2((std::string )"Yb", (uint16_t)2, (double[]){0.0855,18.512,0.2943,7.373,0.6412,2.678,-0.0213,0.0421}, // <j0>
        (double[]){0.3852,18.55,0.9415,6.551,0.8492,2.043,0.0301,0.0478}, // <j2>
        (double[]){-0.2393,9.947,0.0663,4.823,0.2009,0.965,0.0122,0.0311}, // <j4>
        (double[]){-0.0739,5.031,0.014,2.03,0.0351,0.508,0.0174,0.0434}); // <j6>
    ion_map["Yb2"]=Yb2;
    static const MagneticIon Yb3((std::string )"Yb", (uint16_t)3, (double[]){0.0416,16.095,0.2849,7.834,0.6961,2.672,-0.0229,0.0344}, // <j0>
        (double[]){0.157,18.555,0.8484,6.54,0.888,2.037,0.0318,0.0498}, // <j2>
        (double[]){-0.2121,8.197,0.0325,3.153,0.1975,0.884,0.0093,0.0435}, // <j4>
        (double[]){-0.0345,5.007,-0.0677,2.02,0.0985,0.549,-0.0076,0.0359}); // <j6>
    ion_map["Yb3"]=Yb3;
    static const MagneticIon U3((std::string )"U", (uint16_t)3, (double[]){0.5058,23.288,1.3464,7.003,-0.8724,4.868,0.0192,0.1507}, // <j0>
        (double[]){4.1582,16.534,2.4675,5.952,-0.0252,0.765,0.0057,0.0822}, // <j2>
        (double[]){-0.9859,16.601,0.6116,6.515,0.602,2.597,-0.001,0.0599}, // <j4>
        (double[]){-0.3797,9.953,0.0459,5.038,0.2748,1.607,0.0016,0.0345}); // <j6>
    ion_map["U3"]=U3;
    static const MagneticIon U4((std::string )"U", (uint16_t)4, (double[]){0.3291,23.548,1.0836,8.454,-0.434,4.12,0.0214,0.1757}, // <j0>
        (double[]){3.7449,13.894,2.6453,4.863,-0.5218,3.192,0.0009,0.0928}, // <j2>
        (double[]){-1.054,16.605,0.4339,6.512,0.6746,2.599,-0.0011,0.0471}, // <j4>
        (double[]){-0.1793,11.896,-0.2269,5.428,0.3291,1.701,0.003,0.0472}); // <j6>
    ion_map["U4"]=U4;
    static const MagneticIon U5((std::string )"U", (uint16_t)5, (double[]){0.365,19.804,3.2199,6.282,-2.6077,5.301,0.0233,0.175}, // <j0>
        (double[]){3.0724,12.546,2.3076,5.231,-0.0644,1.474,0.0035,0.0477}, // <j2>
        (double[]){-0.9588,16.485,0.1576,6.44,0.7785,2.64,-0.001,0.0493}, // <j4>
        (double[]){-0.0399,11.891,-0.3458,5.58,0.334,1.645,0.0029,0.0444}); // <j6>
    ion_map["U5"]=U5;
    static const MagneticIon Np3((std::string )"Np", (uint16_t)3, (double[]){0.5157,20.865,2.2784,5.893,-1.8163,4.846,0.0211,0.1378}, // <j0>
        (double[]){3.717,15.133,2.3216,5.503,-0.0275,0.8,0.0052,0.0948}, // <j2>
        (double[]){0.9029,16.586,0.4006,6.47,0.6545,2.563,-0.0004,0.047}, // <j4>
        (double[]){-0.2427,11.844,-0.1129,5.377,0.2848,1.568,0.0022,0.0368}); // <j6>
    ion_map["Np3"]=Np3;
    static const MagneticIon Np4((std::string )"Np", (uint16_t)4, (double[]){0.4206,19.805,2.8004,5.978,-2.2436,4.985,0.0228,0.1408}, // <j0>
        (double[]){2.9203,14.646,2.5979,5.559,-0.0301,0.367,0.0141,0.0532}, // <j2>
        (double[]){-0.9887,12.441,0.5918,5.294,0.5306,2.263,-0.0021,0.0583}, // <j4>
        (double[]){-0.2436,9.599,-0.1317,4.101,0.3029,1.545,0.0019,0.05}); // <j6>
    ion_map["Np4"]=Np4;
    static const MagneticIon Np5((std::string )"Np", (uint16_t)5, (double[]){0.3692,18.19,3.151,5.85,-2.5446,4.916,0.0248,0.1515}, // <j0>
        (double[]){2.3308,13.654,2.7219,5.494,-0.1357,0.049,0.1224,0.0553}, // <j2>
        (double[]){-0.8146,16.581,-0.0055,6.475,0.7956,2.562,-0.0004,0.06}, // <j4>
        (double[]){-0.1157,9.565,-0.2654,4.26,0.3298,1.549,0.0025,0.0495}); // <j6>
    ion_map["Np5"]=Np5;
    static const MagneticIon Np6((std::string )"Np", (uint16_t)6, (double[]){0.2929,17.561,3.4866,5.785,-2.8066,4.871,0.0267,0.1698}, // <j0>
        (double[]){1.8245,13.18,2.8508,5.407,-0.1579,0.044,0.1438,0.0585}, // <j2>
        (double[]){0.6738,16.553,-0.2297,6.505,0.8513,2.553,-0.0003,0.0623}, // <j4>
        (double[]){-0.0128,9.569,-0.3611,4.304,0.3419,1.541,0.0032,0.052}); // <j6>
    ion_map["Np6"]=Np6;
    static const MagneticIon Pu3((std::string )"Pu", (uint16_t)3, (double[]){0.384,16.679,3.1049,5.421,-2.5148,4.551,0.0263,0.128}, // <j0>
        (double[]){2.0885,12.871,2.5961,5.19,-0.1465,0.039,0.1343,0.0866}, // <j2>
        (double[]){-0.7014,16.369,-0.1162,6.697,0.7778,2.45,0,0.0546}, // <j4>
        (double[]){-0.0364,9.572,-0.3181,4.342,0.321,1.523,0.0041,0.0496}); // <j6>
    ion_map["Pu3"]=Pu3;
    static const MagneticIon Pu4((std::string )"Pu", (uint16_t)4, (double[]){0.4934,16.836,1.6394,5.638,-1.1581,4.14,0.0248,0.1242}, // <j0>
        (double[]){2.7244,12.926,2.3387,5.163,-0.13,0.046,0.1177,0.049}, // <j2>
        (double[]){-0.916,12.203,0.4891,5.127,0.529,2.149,-0.0022,0.052}, // <j4>
        (double[]){-0.2394,7.837,-0.0785,4.024,0.2643,1.378,0.0012,0.0414}); // <j6>
    ion_map["Pu4"]=Pu4;
    static const MagneticIon Pu5((std::string )"Pu", (uint16_t)5, (double[]){0.3888,16.559,2.0362,5.657,-1.4515,4.255,0.0267,0.1287}, // <j0>
        (double[]){2.1409,12.832,2.5664,5.152,-0.1338,0.046,0.121,0.0491}, // <j2>
        (double[]){-0.7035,16.36,-0.0979,6.706,0.7726,2.447,0,0.061}, // <j4>
        (double[]){-0.109,7.819,-0.2243,4.1,0.2947,1.404,0.0015,0.0477}); // <j6>
    ion_map["Pu5"]=Pu5;
    static const MagneticIon Pu6((std::string )"Pu", (uint16_t)6, (double[]){0.3172,16.051,3.4654,5.351,-2.8102,4.513,0.0281,0.1382}, // <j0>
        (double[]){1.7262,12.324,2.6652,5.066,-0.1695,0.041,0.155,0.0502}, // <j2>
        (double[]){-0.556,16.322,-0.3046,6.768,0.8146,2.426,0.0001,0.0596}, // <j4>
        (double[]){-0.0001,7.82,-0.3354,4.144,0.3097,1.403,0.002,0.0513}); // <j6>
    ion_map["Pu6"]=Pu6;
    static const MagneticIon Am2((std::string )"Am", (uint16_t)2, (double[]){0.4743,21.776,1.58,5.69,-1.0779,4.145,0.0218,0.1253}, // <j0>
        (double[]){3.5237,15.955,2.2855,5.195,-0.0142,0.585,0.0033,0.112}, // <j2>
        (double[]){-0.7433,16.416,0.3481,6.788,0.6014,2.346,0,0.0566}, // <j4>
        (double[]){-0.3176,7.864,0.0771,4.161,0.2194,1.339,0.0018,0.0374}); // <j6>
    ion_map["Am2"]=Am2;
    static const MagneticIon Am3((std::string )"Am", (uint16_t)3, (double[]){0.4239,19.574,1.4573,5.872,-0.9052,3.968,0.0238,0.1054}, // <j0>
        (double[]){2.8622,14.733,2.4099,5.144,-0.1326,0.031,0.1233,0.0727}, // <j2>
        (double[]){0.8092,12.854,0.4161,5.459,0.5476,2.172,-0.0011,0.053}, // <j4>
        (double[]){-0.3159,6.982,0.0682,3.995,0.2141,1.188,-0.0015,0.0281}); // <j6>
    ion_map["Am3"]=Am3;
    static const MagneticIon Am4((std::string )"Am", (uint16_t)4, (double[]){0.3737,17.862,1.3521,6.043,-0.7514,3.72,0.0258,0.1113}, // <j0>
        (double[]){2.4141,12.948,2.3687,4.945,-0.249,0.022,0.2371,0.0502}, // <j2>
        (double[]){-0.8548,12.226,0.3037,5.909,0.6173,2.188,-0.0016,0.0456}, // <j4>
        (double[]){-0.1787,7.88,-0.1274,4.09,0.2565,1.315,0.0017,0.0419}); // <j6>
    ion_map["Am4"]=Am4;
    static const MagneticIon Am5((std::string )"Am", (uint16_t)5, (double[]){0.2956,17.372,1.4525,6.073,-0.7755,3.662,0.0277,0.1202}, // <j0>
        (double[]){2.0109,12.053,2.4155,4.836,-0.2264,0.027,0.2128,0.0414}, // <j2>
        (double[]){-0.6538,15.462,-0.0948,5.997,0.7295,2.297,0,0.0594}, // <j4>
        (double[]){-0.0927,6.073,-0.2227,3.784,0.2916,1.372,0.0026,0.0485}); // <j6>
    ion_map["Am5"]=Am5;
    static const MagneticIon Am6((std::string )"Am", (uint16_t)6, (double[]){0.2302,16.953,1.4864,6.116,-0.7457,3.543,0.0294,0.1323}, // <j0>
        (double[]){1.6778,11.337,2.4531,4.725,-0.2043,0.034,0.1892,0.0387}, // <j2>
        (double[]){-0.539,15.449,-0.2689,6.017,0.7711,2.297,0.0002,0.0729}, // <j4>
        (double[]){0.0152,6.079,-0.3549,3.861,0.3125,1.403,0.0036,0.0732}); // <j6>
    ion_map["Am6"]=Am6;
    static const MagneticIon Am7((std::string )"Am", (uint16_t)7, (double[]){0.3601,12.73,1.964,5.12,-1.356,3.714,0.0316,0.1232}, // <j0>
        (double[]){1.8845,9.161,2.0746,4.042,-0.1318,1.723,0.002,0.0379}, // <j2>
        (double[]){-0.4688,12.019,-0.2692,7.042,0.7297,2.164,-0.0011,0.0262}, // <j4>
        (double[]){0.1292,6.082,-0.4689,3.879,0.3234,1.393,0.0042,0.0475}); // <j6>
    ion_map["Am7"]=Am7;
  }
  return 0;
}

MagneticIon getMagneticIon(const std::string symbol,const uint16_t charge)
{
  initializeMap();
  std::stringstream what;
  what<<symbol<<charge;
  std::map<std::string,MagneticIon>::iterator it;
  it=ion_map.find(what.str());   //check if there is such symbol+charge combination

  if (it==ion_map.end())
  {
    //no such combination 
    std::stringstream msg;
    msg << "Failed to find an atom with symbol=" << symbol << " and charge=" << charge;
    throw std::runtime_error(msg.str());
  }
  else
  {
    return it->second;
  }
}

std::vector <double> getJL(const std::string symbol,const uint16_t charge, const uint16_t l)
{
  MagneticIon temp(getMagneticIon(symbol,charge));
  std::vector <double> v(8,0.);  
  switch (l)
    {
      case 0: {v=temp.j0; break;}
      case 2: {v=temp.j2; break;}
      case 4: {v=temp.j4; break;}
      case 6: {v=temp.j6; break;}
      default: {
          // all other <jl> are not defined    
          std::stringstream msg;
          msg << "Failed to find <j" << l << ">";
          throw std::runtime_error(msg.str());
          }//error
    }
  return v;  
}

} // namespace PhysicalConstants
} // namespace Mantid

