#ifndef _GRE_Util_Noncopyable_HPP
#define _GRE_Util_Noncopyable_HPP

namespace GRE{
	namespace Util{ // protection from unintended ADL
	//  Private copy constructor and copy assignment ensure classes derived from
	//  class noncopyable cannot be copied.
	//  Contributed by Dave Abrahams
  	class Noncopyable
  	{
   	protected:
		Noncopyable() {}
      	~Noncopyable() {}
   	private:  // emphasize the following members are private
      	Noncopyable( const Noncopyable& );
      	const Noncopyable& operator=( const Noncopyable& );
  	};
}//namespace-util
typedef Util::Noncopyable Noncopyable;
}
#endif  
