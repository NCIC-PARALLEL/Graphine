#ifndef GRE_GraphTypesDef_HPP
#define GRE_GraphTypesDef_HPP
namespace GRE {
	typedef uint32_t lvid_t;
	typedef size_t lidx_t;
	typedef uint64_t vid_t;
	typedef size_t idx_t;
	typedef std::pair<lvid_t, lvid_t> ledge_t;
	typedef std::pair<vid_t, vid_t> edge_t;

	#define NotALvid (-1)
}
#endif
