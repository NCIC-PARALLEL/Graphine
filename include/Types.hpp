#ifndef _Types_HPP
#define _Types_HPP
//#include "./Graph/GraphTypesDef.hpp"
extern "C"{
#include <stdint.h>
}
#include <cstddef>
#include <utility>
#include <limits>
//define as visible globally
typedef char byte;

typedef uint32_t lvid_t;
typedef size_t lidx_t;
typedef uint64_t vid_t;
typedef size_t idx_t;

typedef size_t leid_t;
typedef std::pair<lvid_t, lvid_t> ledge_t;
typedef std::pair<vid_t, vid_t> edge_t;

typedef int partID_t;
typedef int subgraph_id_t;

template <typename T>
struct WeightedEdge{
	edge_t edge;
	T property;
};
#define NotALvid 0xffffffffL

enum Dir{Out, In, InOut, None};
namespace GRE{
	//Note this type is not zero-cost.
	typedef struct notUsed{} notUsed_t;
	//define byte as primitive type
}
#endif
