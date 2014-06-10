#ifndef GRE_Message_HPP
#define GRE_Message_HPP
#include "Types.hpp"
namespace GRE{
	template <typename data_t, typename dest_t = vid_t>
	struct Msg{
		dest_t dest;
		data_t data;
		public:
			Msg(dest_t _dest, data_t _data):
				dest(_dest), data(_data){}
			Msg(){}
			~Msg(){}
	};
}
#endif
