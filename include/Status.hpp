#ifndef _GRE_Status_HPP
#define _GRE_Status_HPP
#include <string>
//Define in global namespace
enum GRE_Status {
	GRE_SUCCESS = 0,
	GRE_FAILURE,
	GRE_NOBUFFER,
	GRE_NULL,
	GRE_FINISHED,
	GRE_UNFINISHED,
	GRE_INITIALIZED,
	GRE_UNINITIALIZED,
	GRE_NotImplemented//end: 7
};

class GRE_Status_Parser {
public:
	std::string parse(GRE_Status status){
		static const std::string words[] = {
			"sucess",//0
			"failre",
			"no_buffer",
			"nothing",
			"finished",
			"unfinished",
			"initialized",
			"uninitialized",
			"not_implemented"//end: 7
		};
		//check status is legal!
		return words[status];
	}
};
#endif
