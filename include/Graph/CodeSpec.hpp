#ifndef GRE_Code_Spec_HPP
#define GRE_Code_Spec_HPP
//VID = AgentMask(2 bits)|MacineID(18bits)|Index(44bits)
namespace GRE{
#define MachineIDOffset 44
#define AgentOffset 62
#define IndexMask 0xfffffffffffUL
#define MachineIDMask 0x3ffff00000000000UL
#define ScatterMask		0x4000000000000000UL
#define CombinerMask	0x8000000000000000UL
#define NotAnIndex IndexMask

#define isScatter(v) (ScatterMask & (v))
#define isCombiner(v) (CombinerMask & (v))
#define getAgentFlag(v) ((v) >> AgentOffset)
#define getIndex(v) (IndexMask & (v))
#define setMachineID(v , m) ((v)|((m)<<MachineIDOffset))
#define getMachineID(v) (((v) & MachineIDMask)>>MachineIDOffset)
#define setScatter(v) (ScatterMask | (v))
#define setCombiner(v) (CombinerMask | (v))
}
#endif
