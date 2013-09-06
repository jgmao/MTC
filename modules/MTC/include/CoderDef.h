#ifndef CODER_DEF_H
#define CODER_DEF_H


#ifndef CODING_MODE
#define CODING_MODE
namespace mtc{
enum class CodingMode : int { CODING_MODE_PQI,
						   CODING_MODE_TPQI,
						   CODING_MODE_JPEG,
						   CODING_MODE_MTC,//MTC
						   CODING_MODE_POST_TP,
               CODING_MODE_TPSS //thin plate spline smoothing
};
#endif

//#include "HuffTree.h"
#include "QGrid.h"
#include <fstream>
static const char* month[] = {"Jan","Feb","Mar","Apr","May","Jun","Jul","Aug","Sep","Oct","Nov","Dec"};


//parameters
enum UNIT_TEST { HAS_UNIT_TEST, NO_UNIT_TEST};

enum class PostBlendingType : int { NO_POST_BLENDING, POST_BLENDING_ONLINE, POST_BLENDING_OFFLINE};

enum class LightingCorrectionType : int { NO_LIGHTING_CORRECTION, HAS_LIGHTING_CORRECTION, PQI_LF_ENCODING, PREDEF_LIGHTING, POISSON_LC};
//#define INT_RST

//init_block size
//static const Size3 _init_tp_block_size(32,32,1);
//static const double _init_pqi_q_size = 
#ifndef jpeg_type
#define jpeg_type
enum class JPEGType  : int {		JPEG_BASELINE,
									JPEG_WASTON,
									JPEG_ADAPTIVE,
									JPEG_LF_3,
									JPEG_LF_6,
									JPEG_CUSTOM_TBL};
#endif

//static std::map<int,string> Huffman_DC_Table;
//static std::map<int,string> Huffman_AC_Table;
//static HuffTree Huffman_DC_Tree;
//static HuffTree Huffman_AC_Tree;
#ifndef ADAPTIVE_TPSS
#define ADAPTIVE_TPSS
#endif

}
#endif


