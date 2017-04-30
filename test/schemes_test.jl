
using AcFDTD
using Base.Test

schm = IISO()
show(schm)
schm = SLF()
show(schm)
schm = IWB()
show(schm)

env = FDTDEnv(0.1,IISO()) 
show(env)
env = FDTDEnv(SLF(),1e4) 
show(env)
geo = CuboidRoom(1.,2.,3.,50*ones(6), env)
show(geo)
geo = CuboidRoom(1.,2.,3.,50*ones(6), FDTDEnv(0.1,IISO()))
show(geo)
#geo = LShapedRoom(1.,2.,3.,0.5,1.5,50*ones(6), FDTDEnv(0.1,IISO()))
#show(geo)
#geo = LShapedRoom(10,20,14,5,18,50*ones(6), FDTDEnv(0.1,IISO()))
#show(geo)

size(geo)

