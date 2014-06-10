#include <iostream>
#include <string>
struct IngressOptions{
			//internal ds definitions
			enum Format{bin, ascii};
			enum Method{edgelists, partitions};
			//internal fields
			Format format;
			Method method;
			//
			IngressOptions():format(ascii),method(edgelists){}
			IngressOptions(Format _format, Method _method):format(_format),method(_method){}
			void setFormat(const std::string sformat){
				switch(sformat){
				case "bin": format=bin;break;
				case "ascii": format=ascii;break;
				default:format=ascii;
				}
			}
			void setMethod(std::string smethod){
				if(smethod=="partitions") method=partitions;
				else method=edgelists;
			}
			bool is_binary() const {return format==bin;}
			bool is_ascii() const {return format==ascii;}
			bool is_edgelists() const {return method==edgelists;}
			bool is_partitions() const {return method==partitions;}
};
int main(){
	IngressOptions ops;
	ops.setFormat("bin");
	ops.setMethod("partitions");
	std::cout<<ops.format<<std::endl;
	std::cout<<ops.method<<std::endl;
	if(ops.is_binary()) std::cout<<"binary format"<<std::endl;
	if(ops.is_edgelists()) std::cout<<"edgelists format"<<std::endl;
	return 0;
}
